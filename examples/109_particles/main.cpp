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
    std::vector<bool> is_inner; // particle is inside block
    std::vector<Scal> owner_init; // block owning particles after initialization
    std::vector<Scal> owner_comm; // block owning particles after communication
    std::vector<std::pair<std::string, std::vector<Scal>>> csvdata;
    Scal sum_weight; // averaging weight, number of particles
    Scal sum_dist; // sum of distances from domain center
  } * ctx(sem);
  auto& t = *ctx;

  auto dump = [&](std::string path){
    if (sem()) {
      t.csvdata.clear();
      t.owner_comm.resize(t.x.size(), m.GetId());
      // Add coordinates
      for (auto d : m.dirs) {
        const std::string name(1, M::Dir(d).letter());
        t.csvdata.emplace_back(name, std::vector<Scal>());
        for (size_t i = 0; i < t.x.size(); ++i) {
          if (t.is_inner[i]) {
            t.csvdata.back().second.push_back(t.x[i][d]);
          }
        }
      }
      // Add other fields
      t.csvdata.emplace_back("owner_init", std::vector<Scal>());
      for (size_t i = 0; i < t.x.size(); ++i) {
        if (t.is_inner[i]) {
          t.csvdata.back().second.push_back(t.owner_init[i]);
        }
      }
      t.csvdata.emplace_back("owner_comm", std::vector<Scal>());
      for (size_t i = 0; i < t.x.size(); ++i) {
        if (t.is_inner[i]) {
          t.csvdata.back().second.push_back(t.owner_comm[i]);
        }
      }
    }
    if (sem.Nested()) {
      DumpCsv(t.csvdata, path, m);
    }
  };

  if (sem()) {
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
    t.owner_init.resize(t.x.size(), m.GetId());
    t.velocity.resize(t.x.size());
    t.is_inner.resize(t.x.size(), true);
    m.flags.particles_halo_radius = var.Double["cutoff"] / m.GetCellSize()[0];
  }
  dump("particles_0.csv");
  for (int step = 0; step < var.Int["nsteps"]; ++step) {
    if (sem()) {
      typename M::CommPartRequest req;
      req.x = &t.x;
      req.attr_scal = {&t.owner_init};
      req.attr_vect = {&t.velocity};
      m.CommPart(req);
    }
    if (sem()) {
      t.is_inner.resize(t.x.size());
      for (size_t i = 0; i < t.x.size(); ++i) {
        t.is_inner[i] = m.IsInnerPoint(t.x[i]);
      }
      // Compute gravity force on inner particles within cutoff
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
              if (dist < cutoff * 0.5) {
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
          t.velocity[i] += dt * t.force[i];
          t.x[i] += dt * t.velocity[i];
        }
      }

      t.sum_weight = 0;
      t.sum_dist = 0;
      for (size_t i = 0; i < t.x.size(); ++i) {
        if (t.is_inner[i]) {
          t.sum_weight += 1;
          t.sum_dist += t.x[i].dist(m.GetGlobalLength() * 0.5);
        }
      }
      m.Reduce(&t.sum_weight, Reduction::sum);
      m.Reduce(&t.sum_dist, Reduction::sum);
    }
    if (sem()) {
      if (m.IsRoot()) {
        std::cout << util::Format(
            "step={} dist={}\n", step, t.sum_dist / t.sum_weight);
      }
    }
  }
  dump("particles_1.csv");
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
  parser.AddVariable<int>("--nsteps", 50).Help("Number of time steps to make");
  parser.AddVariable<int>("--ncirc", 64)
      .Help("Number of particles in along circle");
  parser.AddVariable<int>("--naxis", 32).Help(
      "Number of particles along axis in z-direction");
  parser.AddVariable<double>("--cutoff", 0.25).Help("Cutoff radius of force");
  parser.AddVariable<double>("--dt", 0.0005).Help("Time step");
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
