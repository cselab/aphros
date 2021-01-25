// Created by Petr Karnakov on 05.01.2020
// Copyright 2020 ETH Zurich

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
#include <dump/vtk.h>
#include <func/init.h>
#include <solver/curv.h>
#include <solver/reconst.h>
#include <solver/vofm.h>
#include <util/hydro.h>

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
void ReadPlain(std::string path, FieldCell<Scal>& u, GMIdx<3>& size) {
  std::ifstream dat(path);
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
    const Multi<const FieldCell<Scal>*> fck, const MapCondFaceFluid& mff,
    M& m) {
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
        const Scal k = (std::abs(um - 0.5) < std::abs(up - 0.5) ? km : kp);
        const Scal v = (up - um) * k * m.GetArea(f);
        if (!IsNan(v)) {
          ffv[f] += v;
        }
      }
    }
  }
  if (sem.Nested("proj")) {
    ProjectVolumeFlux(ffv, mff, m);
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

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();

  struct {
    std::unique_ptr<Vofm<M>> as; // advection solver
    FieldCell<Scal> fcs; // volume source
    FieldFace<Scal> ffv; // volume flux
    Multi<FieldCell<Scal>> fck; // curvature
    MapCondFaceFluid mf_cond_fluid; // face conditions
    MapCondFaceAdvection<Scal> mf_cond; // face conditions
    GRange<size_t> layers;
    typename PartStrMeshM<M>::Par psm_par;
    std::unique_ptr<PartStrMeshM<M>> psm;
    Scal dt = 0.01;
    size_t step = 0;
    size_t frame = 0;
  } * ctx(sem);
  auto& fcs = ctx->fcs;
  auto& ffv = ctx->ffv;
  auto& fck = ctx->fck;
  auto& layers = ctx->layers;
  auto& psm_par = ctx->psm_par;
  auto& psm = ctx->psm;
  auto& dt = ctx->dt;
  auto& step = ctx->step;
  auto& frame = ctx->frame;

  auto& as = ctx->as;
  const Scal tmax = 10;
  const Scal dtmax = 0.1;
  const Scal cfl = 0.5;

  if (sem("init")) {
    fcs.Reinit(m, 0);
    ffv.Reinit(m, 0);

    GetFluidFaceCond(var, m, ctx->mf_cond_fluid, ctx->mf_cond);

    typename Vofm<M>::Par p;
    p.sharpen = true;
    p.sharpen_cfl = 0.1;
    layers = GRange<size_t>(p.layers);

    Multi<FieldCell<Scal>> fccl; // initial color
    Multi<FieldCell<Scal>> fcu; // initial volume fraction
    const std::string path = "ref/voronoi/points1.dat";
    ReadColorPlain(path, layers, fcu, fccl, m);

    as.reset(new Vofm<M>(m, fcu, fccl, ctx->mf_cond, &ffv, &fcs, 0., dt, p));
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
    psm = UCurv<M>::CalcCurvPart(
        layers, as->GetAlpha(), as->GetNormal(), as->GetMask(), as->GetColor(),
        psm_par, fck, m);
  }
  if (sem.Nested("flux")) {
    CalcMeanCurvatureFlowFlux(
        ffv, layers, as->GetFieldM(), as->GetColor(), as->GetNormal(),
        as->GetAlpha(), fck, ctx->mf_cond_fluid, m);
  }
  if (sem("dt")) {
    Scal maxv = 0;
    for (auto f : m.Faces()) {
      maxv = std::max(maxv, std::abs(ffv[f]));
    }
    dt = cfl * m.GetCellSize().prod() / maxv;
    m.Reduce(&dt, "min");
  }
  if (sem("mindt")) {
    dt = std::min(dt, dtmax);
    as->SetTimeStep(dt);
    if (m.IsRoot()) {
      std::cout << "dt=" << dt << std::endl;
    }
  }
  const size_t dumpskip = 100;
  const bool dump = (step % dumpskip == 0);
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
  const int bsx = 32;
  int nx = 64; // mesh size
  if (argc > 1) {
    nx = atoi(argv[1]);
  }
  int px = 1; // number of ranks in x (assuming py=px)
  if (argc > 2) {
    px = atoi(argv[2]);
  }
  assert(nx % (bsx * px) == 0);
  const int bx = nx / (bsx * px);
  std::string conf = R"EOF(
set int bx 1
set int by 1
set int bz 1

set int bsx 32
set int bsy 32
set int bsz 1

set int dim 2

set int px 8
set int py 8
set int pz 1

set double bcc_clear0 0
set double bcc_clear1 1
set double inletcl 0
set double bcc_fill -1
set string bc_xm symm
set string bc_xp symm
set string bc_ym symm
set string bc_yp symm

set double hypre_symm_tol 1e-10
set int hypre_symm_maxiter 1000
set int hypre_periodic_x 0
set int hypre_periodic_y 0
)EOF";

  conf += "set int bx " + std::to_string(bx) + "\n";
  conf += "set int by " + std::to_string(bx) + "\n";
  conf += "set int px " + std::to_string(px) + "\n";
  conf += "set int py " + std::to_string(px) + "\n";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
