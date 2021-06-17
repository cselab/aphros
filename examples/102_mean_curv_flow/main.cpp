// Created by Petr Karnakov on 28.12.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>

#include <debug/isnan.h>
#include <distr/distrbasic.h>
#include <dump/dumper.h>
#include <dump/vtk.h>
#include <func/init.h>
#include <func/init_u.h>
#include <linear/linear.h>
#include <parse/argparse.h>
#include <solver/curv.h>
#include <solver/reconst.h>
#include <solver/vofm.h>
#include <util/hydro.h>
#include <util/linear.h>
#include "func/init_bc.h"

using M = MeshCartesian<double, 2>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
using R = Reconst<Scal>;
constexpr Scal kClNone = Vofm<M>::kClNone;

// Computes flux proportional to curvature.
// Output:
// ff_flux: result
template <class M>
void CalcMeanCurvatureFlowFlux(
    FieldFace<Scal>& ff_flux, const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*> fcu,
    const Multi<const FieldCell<Scal>*> fccl,
    const Multi<const FieldCell<Vect>*> fcn,
    const Multi<const FieldCell<Scal>*> fca,
    const Multi<const FieldCell<Scal>*> fck,
    const MapEmbed<BCondFluid<Vect>>& mebc, Scal gamma, bool divfree,
    Scal voidpenal, std::shared_ptr<linear::Solver<M>> linsolver, M& m) {
  (void)fcn;
  (void)fca;
  auto sem = m.GetSem();

  if (sem("calc")) {
    ff_flux.Reinit(m, 0);
    for (auto f : m.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      std::set<Scal> s;
      for (auto i : layers) {
        Scal clm = (*fccl[i])[cm];
        Scal clp = (*fccl[i])[cp];
        if (clm != kClNone) s.insert(clm);
        if (clp != kClNone) s.insert(clp);
      }
      for (auto cl : s) {
        Scal um = 0;
        Scal up = 0;
        Scal km = 0;
        Scal kp = 0;
        for (auto i : layers) {
          if ((*fccl[i])[cm] == cl) {
            um = (*fcu[i])[cm];
            km = (*fck[i])[cm];
          }
          if ((*fccl[i])[cp] == cl) {
            up = (*fcu[i])[cp];
            kp = (*fck[i])[cp];
          }
        }
        if (!IsNan(km) || !IsNan(kp)) {
          const Scal k = (IsNan(km) ? kp : IsNan(kp) ? km : (km + kp) * 0.5);
          ff_flux[f] += (up - um) * k * m.GetArea(f) * gamma;
        }
      }
    }
  }
  if (divfree && sem.Nested("proj")) {
    ProjectVolumeFlux(ff_flux, mebc, linsolver, m);
  }
  if (voidpenal && sem("void-penal")) {
    for (auto f : m.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      Scal um = 0;
      Scal up = 0;
      for (auto l : layers) {
        if ((*fccl[l])[cm] != kClNone) {
          um += (*fcu[l])[cm];
        }
        if ((*fccl[l])[cp] != kClNone) {
          up += (*fcu[l])[cp];
        }
      }
      um = std::min(1., um);
      up = std::min(1., up);
      ff_flux[f] += -(up - um) * voidpenal * m.GetArea(f);
    }
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    std::unique_ptr<Vofm<M>> as; // advection solver
    FieldCell<Scal> fc_src; // volume source
    FieldEmbed<Scal> fe_flux; // volume flux
    FieldCell<Vect> fc_vel; // velocity
    Multi<FieldCell<Scal>> fck; // curvature
    MapEmbed<BCondFluid<Vect>> mebc_fluid; // face conditions
    MapEmbed<BCondAdvection<Scal>> mebc_adv; // face conditions
    GRange<size_t> layers;
    typename PartStrMeshM<M>::Par psm_par;
    std::unique_ptr<PartStrMeshM<M>> psm;
    Scal dt = 0;
    size_t step = 0;
    size_t frame = 0;
    std::unique_ptr<Dumper> dumper;
    std::shared_ptr<linear::Solver<M>> linsolver;
  } * ctx(sem);
  auto& layers = ctx->layers;
  auto& t = *ctx;

  auto& as = ctx->as;
  const Scal tmax = var.Double["tmax"];
  const Scal dtmax = var.Double["dtmax"];
  const Scal cfl = var.Double["cfl"];

  if (sem("init")) {
    t.dumper = std::make_unique<Dumper>(var, "dump_");
    t.fc_src.Reinit(m, 0);
    t.fe_flux.Reinit(m, 0);

    t.linsolver = ULinear<M>::MakeLinearSolver(var, "symm", m);

    m.flags.is_periodic[0] = var.Int["hypre_periodic_x"];
    m.flags.is_periodic[1] = var.Int["hypre_periodic_y"];
    m.flags.linreport = var.Int["linreport"];

    {
      auto p = InitBc(var, m, {}, {});
      ctx->mebc_fluid = std::get<0>(p);
      ctx->mebc_adv = std::get<1>(p);
    }

    typename Vofm<M>::Par p;
    p.sharpen = var.Int["sharpen"];
    p.sharpen_cfl = var.Double["sharpen_cfl"];
    p.extrapolate_boundaries = var.Int["extrapolate_boundaries"];
    p.layers = var.Int["vofm_layers"];
    layers = GRange<size_t>(p.layers);

    Multi<FieldCell<Scal>> fccl(layers); // initial color
    Multi<FieldCell<Scal>> fcu(layers); // initial volume fraction
    auto buf = ReadPrimList(var.String["list_path"], m.IsRoot());
    InitOverlappingComponents(buf, fcu, fccl, layers, m);
    as.reset(new Vofm<M>(
        m, m, fcu, fccl, ctx->mebc_adv, &t.fe_flux, &t.fc_src, 0, t.dt, p));

    auto mod = [var](auto& fcu, auto& fccl, auto layers, auto& m) {
      auto buf = ReadPrimList(var.String["list_path"], m.IsRoot());
      InitOverlappingComponents(buf, fcu, fccl, layers, m);
    };
    as->AddModifier(mod);

    t.fck.Reinit(layers, m, 0);
    t.psm_par.dump_fr = 1;
    t.psm_par.dim = 2;
    t.psm_par.ns = 2;
  }
  sem.LoopBegin();
  if (sem.Nested("start")) {
    as->StartStep();
  }
  if (sem.Nested("iter")) {
    as->MakeIteration();
  }
  if (sem.Nested("finish")) {
    as->FinishStep();
  }
  if (sem.Nested("post")) {
    as->PostStep();
  }
  if (sem.Nested()) {
    t.psm = UCurv<M>::CalcCurvPart(as->GetPlic(), t.psm_par, t.fck, m, m);
  }
  if (sem.Nested("flux")) {
    CalcMeanCurvatureFlowFlux(
        t.fe_flux.GetFieldFace(), layers, as->GetFieldM(), as->GetColor(),
        as->GetNormal(), as->GetAlpha(), t.fck, ctx->mebc_fluid,
        var.Double["gamma"], var.Int["divfree"], var.Double["voidpenal"],
        t.linsolver, m);
  }
  if (sem("dt")) {
    Scal maxv = 0;
    for (auto f : m.Faces()) {
      maxv = std::max(maxv, std::abs(t.fe_flux[f]));
    }
    t.dt = cfl * m.GetCellSize().prod() / maxv;
    m.Reduce(&t.dt, "min");
  }
  if (sem("mindt")) {
    t.dt = std::min(t.dt, dtmax);
    as->SetTimeStep(t.dt);
    if (m.IsRoot()) {
      std::cout << "t=" << as->GetTime() << " dt=" << t.dt << std::endl;
    }
  }
  const bool dump = t.dumper->Try(as->GetTime(), as->GetTimeStep());
  if (sem.Nested() && var.Int["dumppoly"] && dump) {
    as->DumpInterface(GetDumpName("s", ".vtk", t.frame));
  }
  if (sem.Nested() && var.Int["dumppolymarch"] && dump) {
    as->DumpInterfaceMarch(GetDumpName("sm", ".vtk", t.frame));
  }
  if (sem.Nested() && var.Int["dumppart"] && dump) {
    t.psm->DumpParticles(
        as->GetAlpha(), as->GetNormal(), t.frame, as->GetTime());
  }
  if (sem("checkloop")) {
    if (dump) {
      ++t.frame;

      t.fc_vel.Reinit(m, Vect(0));
      for (auto c : m.AllCells()) {
        for (auto q : m.Nci(c)) {
          const auto f = m.GetFace(c, q);
          t.fc_vel[c] += m.GetNormal(f) * (t.fe_flux[f] * 0.5 / m.GetArea(f));
        }
      }
      m.Dump(&t.fc_vel, 0, "vx");
      m.Dump(&t.fc_vel, 1, "vy");
    }
    if (as->GetTime() >= tmax) {
      sem.LoopBreak();
    }
    ++t.step;
  }
  sem.LoopEnd();
}

int main(int argc, const char** argv) {
  ArgumentParser parser("Constrained mean curvature flow");
  parser.AddVariable<std::string>("config", "std.conf")
      .Help("Path to config file");
  parser.AddVariable<int>("--nx", 16).Help("Mesh size in x-direction");
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  parser.AddVariable<int>("--ny", 0).Help(
      "Mesh size in y-direction. Defaults to NX");
  parser.AddVariable<int>("--bs", 16).Help("Block size in all directions");
  const auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  const auto confpath = args.String["config"];
  std::ifstream fconf(confpath);

  if (!fconf.good()) {
    throw std::runtime_error("Can't open config '" + confpath + "'");
  }

  std::stringstream conf;
  conf << fconf.rdbuf();

  const int nx = args.Int["nx"];
  const int ny = args.Int["ny"];
  const MIdx mesh_size(nx, ny ? ny : nx);
  const MIdx block_size(args.Int["bs"], args.Int["bs"]);

  MpiWrapper mpi(&argc, &argv);
  Subdomains<MIdx> sub(mesh_size, block_size, mpi.GetCommSize());
  conf << sub.GetConfig() << '\n';
  conf << args.String["extra"] << '\n';

  return RunMpiBasicString<M>(mpi, Run, conf.str());
}
