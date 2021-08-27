// Created by Petr Karnakov on 16.08.2021
// Copyright 2021 ETH Zurich

#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <distr/distrbasic.h>
#include <dump/hdf.h>
#include <func/init.h>
#include <parse/argparse.h>
#include <parse/vars.h>
#include <util/distr.h>
#include <util/filesystem.h>
#include <util/format.h>
#include <util/vof.h>

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

std::ofstream fout;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::vector<Vect> x;
    std::vector<Vect> velocity;
    std::vector<Vect> force;
    std::vector<bool> is_inner; // particle is owned by current block
    std::vector<Scal> block_init; // block owning particles after initialization
    std::vector<Scal> block_comm; // block owning particles after communication
    std::vector<Scal> rank_init; // rank owning particles after communication
    std::vector<std::pair<std::string, std::vector<Scal>>> csvdata;
    Scal npart; // total number of particles inner particles
    Scal dist_mean; // mean distance from axis
    Scal dist_var; // variance of distance from z-axis
    int step = 0;
  } * ctx(sem);
  auto& t = *ctx;

  auto dump = [&](std::string path, bool only_inner){
    if (sem("dump_local")) {
      t.csvdata.clear();
      t.block_comm.resize(t.x.size(), m.GetId());
      auto append = [&](std::string name, auto func) {
        t.csvdata.emplace_back(name, std::vector<Scal>());
        for (size_t i = 0; i < t.x.size(); ++i) {
          if (!only_inner || t.is_inner[i]) {
            func(t.csvdata.back().second, i);
          }
        }
      };
      // Add positions
      for (auto d : m.dirs) {
        const std::string name(1, M::Dir(d).letter());
        append(name, [&](auto& data, size_t i) { //
          data.push_back(t.x[i][d]);
        });
      }
      // Add other fields
      append("block_init", [&](auto& data, size_t i) { //
        data.push_back(t.block_init[i]);
      });
      append("block_comm", [&](auto& data, size_t i) { //
        data.push_back(t.block_comm[i]);
      });
      append("rank_init", [&](auto& data, size_t i) { //
        data.push_back(m.GetMpiRankFromId(t.block_init[i]));
      });
      append("rank_comm", [&](auto& data, size_t i) { //
        data.push_back(m.GetMpiRankFromId(t.block_comm[i]));
      });
      append("inner", [&](auto& data, size_t i) { //
        data.push_back(t.is_inner[i]);
      });
      if (m.IsRoot()) {
        std::cerr << path << '\n';
      }
    }
    if (sem.Nested()) {
      DumpCsv(t.csvdata, path, m);
    }
  };

  if (sem("init")) {
    const int ncirc = var.Int["ncirc"];
    const int naxis = var.Int["naxis"];
    // Seed particles on cylinder
    for (int ic = 0; ic < ncirc; ++ic) {
      for (int ia = 0; ia < naxis; ++ia) {
        Vect x(0);
        const Scal a = 2 * M_PI * ic / ncirc;
        x[0] = 0.5 + std::sin(a) * 0.25;
        x[1] = 0.5 + std::cos(a) * 0.25;
        x[2] = Scal(ia + 0.5) / naxis;
        if (m.IsInnerPoint(x)) {
          t.x.push_back(x);
        }
      }
    }
    t.block_init.resize(t.x.size(), m.GetId());
    t.velocity.resize(t.x.size());
    t.is_inner.resize(t.x.size(), true);
    m.flags.particles_halo_radius = var.Double["cutoff"] / m.GetCellSize()[0];
  }
  sem.LoopBegin();
  if (sem("commpart")) {
    typename M::CommPartRequest req;
    req.x = &t.x;
    req.is_inner = &t.is_inner;
    req.attr_scal = {&t.block_init};
    req.attr_vect = {&t.velocity};
    m.CommPart(req);
  }
  dump(util::Format("particles_{:04d}.csv", t.step), false);
  if (sem("stat_mean")) {
    t.npart = 0;
    t.dist_mean = 0;
    for (size_t i = 0; i < t.x.size(); ++i) {
      if (t.is_inner[i]) {
        t.npart += 1;
        Vect dx = t.x[i] - m.GetGlobalLength() * 0.5;
        dx[2] = 0;
        t.dist_mean += dx.norm();
      }
    }
    m.Reduce(&t.npart, Reduction::sum);
    m.Reduce(&t.dist_mean, Reduction::sum);
  }
  if (sem("stat_var")) {
    t.dist_mean /= t.npart;
    t.dist_var = 0;
    for (size_t i = 0; i < t.x.size(); ++i) {
      if (t.is_inner[i]) {
        Vect dx = t.x[i] - m.GetGlobalLength() * 0.5;
        dx[2] = 0;
        t.dist_var += sqr(dx.norm() - t.dist_mean);
      }
    }
    m.Reduce(&t.dist_var, Reduction::sum);
  }
  if (sem("stat_print")) {
    t.dist_var /= t.npart;
    if (m.IsRoot()) {
      std::cout << util::Format(
          "step={:02d} npart={} mean={:.8e} std={:.8e}\n", t.step, t.npart,
          t.dist_mean, std::sqrt(t.dist_var));
    }
    if (++t.step > var.Int["nsteps"]) {
      sem.LoopBreak();
    }
  }
  if (sem("advance")) {
    // Compute gravity force for inner particles within cutoff
    const Scal cutoff = var.Double["cutoff"];
    t.force.resize(t.x.size());
    for (size_t i = 0; i < t.x.size(); ++i) {
      if (t.is_inner[i]) {
        auto& f = t.force[i];
        f = Vect(0);
        for (size_t j = 0; j < t.x.size(); ++j) {
          const auto xi = t.x[i];
          const auto xj = t.x[j];
          if (xi != xj) {
            const Scal dist = xi.dist(xj);;
            if (dist < cutoff) {
              f += (xi - xj) / std::pow(dist + 0.01, 3);
            }
          }
        }
      }
    }
    // Apply acceleration and displacement to inner particles
    const Scal dt = var.Double["dt"];
    for (size_t i = 0; i < t.x.size(); ++i) {
      if (t.is_inner[i]) {
        //t.velocity[i] += dt * t.force[i];
        //t.x[i] += dt * t.velocity[i];
        t.x[i] += dt * Vect(3, 2, 1);
      }
    }
  }
  sem.LoopEnd();
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);
  fout.open(util::Format("out_{}", mpi.GetCommRank()));

  ArgumentParser parser(
      "Particles with gravity initialized on a cylinder", mpi.IsRoot());
  auto instances = []() {
    auto map = ModuleDistr<M>::GetInstances();
    std::vector<std::string> res;
    for (auto p : map) {
      res.push_back(p.first);
    }
    return res;
  }();
  parser.AddVariable<std::string>("--backend", "local")
      .Help("Communication backend")
      .Options(instances);
  parser.AddVariable<int>("--block", 16).Help("Block size in all directions");
  parser.AddVariable<int>("--mesh", 32).Help("Mesh size in all directions");
  parser.AddVariable<int>("--nsteps", 10).Help("Number of time steps to make");
  parser.AddVariable<int>("--ncirc", 64)
      .Help("Number of particles in along circle");
  parser.AddVariable<int>("--naxis", 32).Help(
      "Number of particles along axis in z-direction");
  parser.AddVariable<double>("--cutoff", 0.15).Help("Cutoff radius of force");
  parser.AddVariable<double>("--dt", 0.001).Help("Time step");
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::string conf;

  MIdx mesh_size(args.Int["mesh"]);
  MIdx block_size(args.Int["block"]);

  Subdomains<MIdx> sub(mesh_size, block_size, mpi.GetCommSize());
  conf += sub.GetConfig();

  conf += "\nset string backend " + args.String["backend"];
  conf += "\nset int ncirc " + args.Int.GetStr("ncirc");
  conf += "\nset int naxis " + args.Int.GetStr("naxis");
  conf += "\nset int nsteps " + args.Int.GetStr("nsteps");
  conf += "\nset double cutoff " + args.Double.GetStr("cutoff");
  conf += "\nset double dt " + args.Double.GetStr("dt");
  conf += "\n" + args.String["extra"];

  return RunMpiBasicString<M>(mpi, Run, conf);
}
