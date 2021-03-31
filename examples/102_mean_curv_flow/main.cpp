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
#include <linear/linear.h>
#include <parse/argparse.h>
#include <solver/curv.h>
#include <solver/reconst.h>
#include <solver/vofm.h>
#include <util/hydro.h>
#include <util/linear.h>
#include "func/init_bc.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
using R = Reconst<Scal>;
constexpr Scal kClNone = -1;

// Reads data in plain format.
// u: scalar field defined on b
// ndc: index cells
// bc: block cells
// op: output path
// Format:
// <Nx> <Ny> <Nz>
// <data:x=0,y=0,z=0> <data:x=1,y=0,z=0> ...
template <class Scal>
void ReadPlain(std::string path, FieldCell<Scal>& u, generic::MIdx<3>& size) {
  std::ifstream dat(path);
  if (!dat.good()) {
    throw std::runtime_error("ReadPlain: Can't open data file '" + path + "'");
  }
  dat >> size[0] >> size[1] >> size[2];
  GIndex<IdxCell, 3> bc(size);
  u.Reinit(bc, 0);
  for (auto c : bc.Range()) {
    dat >> u[c];
  }
}

template <class M>
void CalcMeanCurvatureFlowFlux(
    FieldFace<Scal>& ffv, const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*> fcu,
    const Multi<const FieldCell<Scal>*> fccl,
    const Multi<const FieldCell<Vect>*> fcn,
    const Multi<const FieldCell<Scal>*> fca,
    const Multi<const FieldCell<Scal>*> fck,
    const MapEmbed<BCondFluid<Vect>>& mff, bool divfree, Scal* voidpenal,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m) {
  (void)fcn;
  (void)fca;
  auto sem = m.GetSem();

  if (sem("calc")) {
    ffv.Reinit(m, 0);
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
          ffv[f] += (up - um) * k * m.GetArea(f);
        }
      }
    }
  }
  if (divfree && sem.Nested("proj")) {
    ProjectVolumeFlux(ffv, mff, linsolver, m);
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
      ffv[f] += -(up - um) * (*voidpenal) * m.GetArea(f);
    }
  }
}

// Returns point at which interpolant has value 0.
// x0,x1: points
// f0,f1: values
static Vect GetIso(Vect x0, Vect x1, Scal f0, Scal f1) {
  return (x0 * f1 - x1 * f0) / (f1 - f0);
}

// cl0: color to filter
// nr: normal
std::vector<Vect> GetPoly(
    const FieldNode<Scal>& fncl, const Scal cl0, const IdxFace f, const M& m) {
  const size_t em = m.GetNumNodes(f);
  std::vector<Vect> xx;
  for (size_t e = 0; e < em; ++e) {
    const size_t ep = (e + 1) % em;
    const IdxNode n = m.GetNode(f, e);
    const IdxNode np = m.GetNode(f, ep);
    const Scal l = (fncl[n] == cl0 ? 1 : -1);
    const Scal lp = (fncl[np] == cl0 ? 1 : -1);
    const Vect x = m.GetNode(n);
    const Vect xp = m.GetNode(np);
    if (l > 0) {
      xx.push_back(x);
    }
    if ((l < 0) != (lp < 0)) {
      xx.push_back(GetIso(x, xp, l, lp));
    }
  }
  return xx;
}

// cl0: color to filter
// nr: normal
Vect GetNormal(
    const FieldNode<Scal>& fncl, const Scal cl0, const IdxCell c, const M& m) {
  Vect n(0);
  for (auto q : m.Nci(c)) {
    const IdxFace f = m.GetFace(c, q);
    const auto xx = GetPoly(fncl, cl0, f, m);
    const Scal area = std::abs(R::GetArea(xx, m.GetNormal(f)));
    n += m.GetNormal(f) * area * m.GetOutwardFactor(c, q);
  }
  return -n / n.norm();
}

// cl0: color to filter
// nr: normal
Scal GetAlpha(
    const FieldNode<Scal>& fncl, const Scal cl0, const IdxCell c, const Vect nr,
    const M& m) {
  Scal a = 0;
  Scal aw = 0;
  for (auto q : m.Nci(c)) {
    const IdxFace f = m.GetFace(c, q);
    const size_t em = m.GetNumNodes(f);
    for (size_t e = 0; e < em; ++e) {
      const size_t ep = (e + 1) % em;
      const IdxNode n = m.GetNode(f, e);
      const IdxNode np = m.GetNode(f, ep);
      const Scal l = (fncl[n] == cl0 ? 1 : -1);
      const Scal lp = (fncl[np] == cl0 ? 1 : -1);
      const Vect x = m.GetNode(n);
      const Vect xp = m.GetNode(np);

      if ((l < 0) != (lp < 0)) {
        a += nr.dot(GetIso(x, xp, l, lp) - m.GetCenter(c));
        aw += 1.;
      }
    }
  }
  return a / aw;
}

template <class M>
void DumpFaces(
    const GRange<size_t>& layers, FieldNode<Scal>& fncl, const M& m) {
  if (m.IsRoot()) {
    std::vector<std::vector<Vect>> vv;
    for (auto c : m.Cells()) {
      std::set<Scal> set;
      for (size_t e = 0; e < m.GetNumNodes(c); ++e) {
        const IdxNode n = m.GetNode(c, e);
        const Scal cl = fncl[n];
        if (cl != kClNone) {
          set.insert(cl);
        }
      }
      for (auto cl : set) {
        for (auto q : m.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          const auto xx = GetPoly(fncl, cl, f, m);
          vv.push_back(xx);
        }
      }
    }
    const std::string fn = "faces.vtk";
    WriteVtkPoly<Vect>(fn, vv, nullptr, {}, {}, "", true, false, false);
  }
}

template <class M>
void InitColorFromNodes(
    const FieldNode<Scal>& fncl, const GRange<size_t>& layers,
    Multi<FieldCell<Scal>>& fcu, Multi<FieldCell<Scal>>& fccl, const M& m) {
  fcu.Reinit(layers, m, 0);
  fccl.Reinit(layers, m, kClNone);

  for (auto c : m.Cells()) {
    std::set<Scal> set;
    // gather colors from adjacent nodes
    for (size_t e = 0; e < m.GetNumNodes(c); ++e) {
      const IdxNode n = m.GetNode(c, e);
      const Scal cl = fncl[n];
      if (cl != kClNone) {
        set.insert(cl);
      }
    }
    // compute volume fraction for every color
    size_t l = 0;
    for (auto cl : set) {
      const Vect nr = GetNormal(fncl, cl, c, m);
      const Scal a = GetAlpha(fncl, cl, c, nr, m);
      if (set.size() == 1) {
        fcu[l][c] = 1;
      } else {
        fcu[l][c] = R::GetLineU(nr, a, m.GetCellSize());
      }
      fccl[l][c] = cl;
      ++l;
      if (l >= layers.size()) {
        break;
      }
    }
  }
}

template <class M>
void ReadColorPlain(
    const std::string path, const GRange<size_t>& layers,
    Multi<FieldCell<Scal>>& fcu, Multi<FieldCell<Scal>>& fccl, const M& m) {
  fcu.resize(layers);
  fccl.resize(layers);
  fcu.InitAll(FieldCell<Scal>(m, 0));
  fccl.InitAll(FieldCell<Scal>(m, kClNone));

  FieldNode<Scal> fncl(m, kClNone);
  FieldCell<Scal> qfccl;
  MIdx qsize;
  ReadPlain(path, qfccl, qsize);
  GIndex<IdxCell, M::dim> qbc(qsize);
  std::cout << "qsize=" << qbc.GetSize() << std::endl;
  auto& bn = m.GetIndexNodes();
  const MIdx size = m.GetGlobalSize() + MIdx(1);
  for (auto n : m.Nodes()) {
    const MIdx w = bn.GetMIdx(n);
    const MIdx qw = w * qsize / size;
    const IdxCell qc = qbc.GetIdx(qw);
    fncl[n] = qfccl[qc];
  }

  InitColorFromNodes(fncl, layers, fcu, fccl, m);
}

template <class M>
void InitColorJunctionT(
    Vect center, const GRange<size_t>& layers, Multi<FieldCell<Scal>>& fcu,
    Multi<FieldCell<Scal>>& fccl, const M& m) {
  FieldNode<Scal> fncl(m, kClNone);
  for (auto n : m.Nodes()) {
    const Vect x = m.GetNode(n);
    const Vect dx = x - center;
    fncl[n] = (dx[1] < 0 ? 0 : dx[0] < 0 ? 1 : 2);
  }

  InitColorFromNodes(fncl, layers, fcu, fccl, m);
}

template <class M>
void InitColorJunctionTSymm(
    Vect center, const GRange<size_t>& layers, Multi<FieldCell<Scal>>& fcu,
    Multi<FieldCell<Scal>>& fccl, const M& m) {
  FieldNode<Scal> fncl(m, kClNone);
  for (auto n : m.Nodes()) {
    const Vect x = m.GetNode(n);
    const Vect dx = x - center;
    const Scal q = x[1] - (1 - center[1]);
    fncl[n] = (dx[1] < 0 ? 0 : q > 0 ? 3 : dx[0] < 0 ? 1 : 2);
  }

  InitColorFromNodes(fncl, layers, fcu, fccl, m);
}

template <class M>
void SetZeroBoundaryFlux(FieldFace<Scal>& fcu, const M& m) {
  for (auto fb : m.Faces()) {
    size_t nci;
    if (m.IsBoundary(fb, nci)) {
      const IdxCell c = m.GetCell(fb, nci);
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        fcu[f] = 0;
      }
    }
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();

  struct {
    std::unique_ptr<Vofm<M>> as; // advection solver
    FieldCell<Scal> fcs; // volume source
    FieldEmbed<Scal> fev; // volume flux
    Multi<FieldCell<Scal>> fck; // curvature
    MapEmbed<BCondFluid<Vect>> mebc_fluid; // face conditions
    MapEmbed<BCondAdvection<Scal>> mebc_adv; // face conditions
    GRange<size_t> layers;
    typename PartStrMeshM<M>::Par psm_par;
    std::unique_ptr<PartStrMeshM<M>> psm;
    Scal dt = 0.01;
    size_t step = 0;
    size_t frame = 0;
    std::unique_ptr<Dumper> dumper;
    std::shared_ptr<linear::Solver<M>> linsolver;
  } * ctx(sem);
  auto& fcs = ctx->fcs;
  auto& fev = ctx->fev;
  auto& fck = ctx->fck;
  auto& layers = ctx->layers;
  auto& psm_par = ctx->psm_par;
  auto& psm = ctx->psm;
  auto& dt = ctx->dt;
  auto& step = ctx->step;
  auto& frame = ctx->frame;
  auto& dumper = ctx->dumper;
  auto& t = *ctx;

  auto& as = ctx->as;
  const Scal tmax = var.Double["tmax"];
  const Scal dtmax = var.Double["dtmax"];
  const Scal cfl = var.Double["cfl"];

  if (sem("init")) {
    dumper = std::make_unique<Dumper>(var, "dump_");
    fcs.Reinit(m, 0);
    fev.Reinit(m, 0);

    t.linsolver = ULinear<M>::MakeLinearSolver(var, "symm", m);

    m.flags.is_periodic[0] = var.Int["hypre_periodic_x"];
    m.flags.is_periodic[1] = var.Int["hypre_periodic_y"];
    m.flags.is_periodic[2] = var.Int["hypre_periodic_z"];

    {
      auto p = InitBc(var, m, {}, {});
      ctx->mebc_fluid = std::get<0>(p);
      ctx->mebc_adv = std::get<1>(p);
    }

    typename Vofm<M>::Par p;
    p.sharpen = var.Int["sharpen"];
    p.sharpen_cfl = var.Double["sharpen_cfl"];
    p.avgnorm0 = var.Double["avgnorm0"];
    p.avgnorm1 = var.Double["avgnorm1"];
    p.extrapolate_boundaries = var.Int["extrapolate_boundaries"];
    layers = GRange<size_t>(p.layers);

    Multi<FieldCell<Scal>> fccl; // initial color
    Multi<FieldCell<Scal>> fcu; // initial volume fraction

    const std::string init_color = var.String["init_color"];
    if (init_color == "triple") {
      const Vect center(var.Vect["triple_center"]);
      InitColorJunctionT(center, layers, fcu, fccl, m);
    } else if (init_color == "triple_symm") {
      const Vect center(var.Vect["triple_center"]);
      InitColorJunctionTSymm(center, layers, fcu, fccl, m);
    } else {
      ReadColorPlain(init_color, layers, fcu, fccl, m);
    }

    as.reset(
        new Vofm<M>(m, m, fcu, fccl, ctx->mebc_adv, &fev, &fcs, 0., dt, p));
    fck.resize(layers);
    fck.InitAll(FieldCell<Scal>(m, 1));
    psm_par.dump_fr = 1;
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
  if (sem.Nested()) {
    psm = UCurv<M>::CalcCurvPart(as->GetPlic(), psm_par, fck, m, m);
  }
  if (sem.Nested("flux")) {
    CalcMeanCurvatureFlowFlux(
        fev.GetFieldFace(), layers, as->GetFieldM(), as->GetColor(),
        as->GetNormal(), as->GetAlpha(), fck, ctx->mebc_fluid,
        var.Int["divfree"], var.Double.Find("voidpenal"), t.linsolver, m);
  }
  if (sem("dt")) {
    if (var.Int["zero_boundary_flux"]) {
      SetZeroBoundaryFlux(fev.GetFieldFace(), m);
    }
    Scal maxv = 0;
    for (auto f : m.Faces()) {
      maxv = std::max(maxv, std::abs(fev[f]));
    }
    dt = cfl * m.GetCellSize().prod() / maxv;
    m.Reduce(&dt, "min");
  }
  if (sem("mindt")) {
    dt = std::min(dt, dtmax);
    as->SetTimeStep(dt);
    if (m.IsRoot()) {
      std::cout << "t=" << as->GetTime() << " dt=" << dt << std::endl;
    }
  }
  const bool dump = dumper->Try(as->GetTime(), as->GetTimeStep());
  if (sem.Nested()) {
    if (dump) {
      as->DumpInterface(GetDumpName("s", ".vtk", frame));
    }
  }
  if (sem.Nested()) {
    if (dump) {
      as->DumpInterfaceMarch(GetDumpName("sm", ".vtk", frame));
    }
  }
  if (sem("checkloop")) {
    if (dump) {
      ++frame;
    }
    if (as->GetTime() >= tmax) {
      sem.LoopBreak();
    }
    ++step;
  }
  sem.LoopEnd();
}

int main(int argc, const char** argv) {
  ArgumentParser parser("Constrained mean curvature flow");
  parser.AddVariable<std::string>("path", "a.conf").Help("Path to config file");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  auto path = args.String["path"];
  std::ifstream fconf(path);

  if (!fconf.good()) {
    throw std::runtime_error("Can't open config '" + path + "'");
  }

  std::stringstream conf;
  conf << fconf.rdbuf();

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf.str());
}
