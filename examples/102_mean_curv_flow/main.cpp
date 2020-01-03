//#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>

#include <debug/isnan.h>
#include <distr/distrbasic.h>
#include <func/init.h>
#include <solver/curv.h>
#include <solver/vofm.h>
#include <solver/reconst.h>
#include <util/hydro.h>

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
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
std::map<Scal, Scal> CalcArea(
    const GRange<size_t>& layers, const Multi<const FieldCell<Vect>*> fcn,
    const Multi<const FieldCell<Scal>*> fca,
    const Multi<const FieldCell<Scal>*> fccl, M& m) {
  auto sem = m.GetSem();
  struct {
    std::vector<Scal> vcl; // color
    std::vector<Scal> vs; // area
  } * ctx(sem);
  auto& vcl = ctx->vcl;
  auto& vs = ctx->vs;
  if (sem("area")) {
    std::map<Scal, Scal> map;
    for (auto c : m.Cells()) {
      for (auto l : layers) {
        const Scal cl = (*fccl[l])[c];
        const Scal a = (*fca[l])[c];
        const Vect n = (*fcn[l])[c];
        if (cl != kClNone && !IsNan(a)) {
          using R = Reconst<Scal>;
          auto xx = R::GetCutPoly2(n, a, m.GetCellSize());
          map[cl] += std::abs(R::GetArea(xx, n));
        }
      }
    }
    for (auto p : map) {
      vcl.push_back(p.first);
      vs.push_back(p.second);
    }
    using T = typename M::template OpCatT<Scal>;
    m.Reduce(std::make_shared<T>(&vcl));
    m.Reduce(std::make_shared<T>(&vs));
  }
  if (sem("bcast")) {
    if (m.IsRoot()) {
      std::map<Scal, Scal> map;
      for (size_t i = 0; i < vcl.size(); ++i) {
        map[vcl[i]] += vs[i];
      }
      vcl.clear();
      vs.clear();
      for (auto p : map) {
        vcl.push_back(p.first);
        vs.push_back(p.second);
      }
    }
    using T = typename M::template OpCatT<Scal>;
    m.Bcast(std::make_shared<T>(&vcl));
    m.Bcast(std::make_shared<T>(&vs));
  }
  if (sem("map")) {
    std::map<Scal, Scal> map;
    for (size_t i = 0; i < vcl.size(); ++i) {
      map[vcl[i]] = vs[i];
    }
    return map;
  }
  return {};
}

template <class M>
std::map<Scal, Scal> CalcVolume(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*> fcu,
    const Multi<const FieldCell<Scal>*> fccl, M& m) {
  auto sem = m.GetSem();
  struct {
    std::vector<Scal> vcl; // color
    std::vector<Scal> vvol; // volume
  } * ctx(sem);
  auto& vcl = ctx->vcl;
  auto& vvol = ctx->vvol;
  if (sem("volume")) {
    std::map<Scal, Scal> map;
    for (auto c : m.Cells()) {
      for (auto l : layers) {
        const Scal cl = (*fccl[l])[c];
        if (cl != kClNone) {
          map[cl] += (*fcu[l])[c] * m.GetVolume(c);
        }
      }
    }
    for (auto p : map) {
      vcl.push_back(p.first);
      vvol.push_back(p.second);
    }
    using T = typename M::template OpCatT<Scal>;
    m.Reduce(std::make_shared<T>(&vcl));
    m.Reduce(std::make_shared<T>(&vvol));
  }
  if (sem("bcast")) {
    if (m.IsRoot()) {
      std::map<Scal, Scal> map;
      for (size_t i = 0; i < vcl.size(); ++i) {
        map[vcl[i]] += vvol[i];
      }
      vcl.clear();
      vvol.clear();
      for (auto p : map) {
        vcl.push_back(p.first);
        vvol.push_back(p.second);
      }
    }
    using T = typename M::template OpCatT<Scal>;
    m.Bcast(std::make_shared<T>(&vcl));
    m.Bcast(std::make_shared<T>(&vvol));
  }
  if (sem("map")) {
    std::map<Scal, Scal> map;
    for (size_t i = 0; i < vcl.size(); ++i) {
      map[vcl[i]] = vvol[i];
    }
    return map;
  }
  return {};
}

template <class M>
void CalcMeanCurvatureFlowFlux(
    FieldFace<Scal>& ffv, const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*> fcu,
    const Multi<const FieldCell<Scal>*> fccl,
    const Multi<const FieldCell<Vect>*> fcn,
    const Multi<const FieldCell<Scal>*> fca,
    const Multi<const FieldCell<Scal>*> fck, const std::map<Scal, Scal>& mvol0,
    const MapCondFaceFluid& mff, M& m) {
  auto sem = m.GetSem();
  struct {
    std::map<Scal, Scal> mvol;
    std::map<Scal, Scal> marea;
  } * ctx(sem);
  auto& mvol = ctx->mvol;
  auto& marea = ctx->marea;

  const Scal kvol = 0.1;
  const Scal kvoid = 0.001;

  if (sem.Nested("volume")) {
    mvol = CalcVolume(layers, fcu, fccl, m);
  }
  if (sem.Nested("area")) {
    marea = CalcArea(layers, fcn, fca, fccl, m);
  }
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
          ffv[f] = v;
        }
      }
    }
  }
  if (sem.Nested("proj")) {
    ProjectVolumeFlux(ffv, mff, m);
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
    std::map<Scal, Scal> mvol0;
  } * ctx(sem);
  auto& fcs = ctx->fcs;
  auto& ffv = ctx->ffv;
  auto& fck = ctx->fck;
  auto& layers = ctx->layers;
  auto& psm_par = ctx->psm_par;
  auto& psm = ctx->psm;
  auto& dt = ctx->dt;
  auto& step = ctx->step;
  auto& mvol0 = ctx->mvol0;

  auto& as = ctx->as;
  const Scal tmax = 10;
  const Scal dtmax = 0.1;
  const Scal cfl = 0.5;

  // if (sem.Nested()) {
  //  InitVf(fcu, var, m);
  //}
  if (sem("init")) {
    fcs.Reinit(m, 0);
    ffv.Reinit(m, 0);

    FieldCell<Scal> fccl(m, kClNone); // initial color
    FieldCell<Scal> fcu(m, 0); // initial volume fraction

    /*
    const std::string qpath = "ref/voronoi/cl.dat";
    FieldCell<Scal> qfccl;
    MIdx qsize;
    ReadPlain(qpath, qfccl, qsize);
    GIndex<IdxCell, M::dim> qbc(qsize);
    std::cout << "qsize=" << qbc.GetSize() << std::endl;
    auto& bc = m.GetIndexCells();
    const MIdx size = m.GetGlobalSize();
    for (auto c : m.Cells()) {
      const MIdx w = bc.GetMIdx(c);
      const MIdx qw = w * qsize / size;
      const IdxCell qc = qbc.GetIdx(qw);
      fccl[c] = qfccl[qc];
    }
    */
    for (auto c : m.AllCells()) {
      auto D = [](Vect x) { return Vect(x[0], x[1], 0);};
      auto rot = [](Vect x) { return Vect(x[0]-x[1]*0.05, x[1]+x[0]*0.05, 0);};
      Vect x = D(m.GetCenter(c));
      x = rot(x);
      fccl[c] = x[1] < 0.3 ? 1 : x[0] < 0.5 ? 2 : 3;
      fcu[c] = 1;
    }

    GetFluidFaceCond(var, m, ctx->mf_cond_fluid, ctx->mf_cond);

    typename Vofm<M>::Par p;
    p.sharpen = true;
    p.sharpen_cfl = 0.1;
    as.reset(new Vofm<M>(m, fcu, fccl, ctx->mf_cond, &ffv, &fcs, 0., dt, p));
    layers = GRange<size_t>(as->GetNumLayers());
    fck.resize(layers);
    fck.InitAll(FieldCell<Scal>(m, 1));
    psm_par.dump_fr = 1;
  }
  if (sem.Nested("volume")) {
    mvol0 = CalcVolume(layers, as->GetFieldM(), as->GetColor(), m);
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
        as->GetAlpha(), fck, mvol0, ctx->mf_cond_fluid, m);
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

    if (step > 1) {
      auto p = as->GetPar();
      // p.sharpen = false;
      as->SetPar(p);
    }
  }
  /*
  if (sem.Nested()) {
    psm->DumpParticles(as->GetAlpha(), as->GetNormal(), step, as->GetTime());
  }
  */
  if (sem.Nested()) {
    as->DumpInterface(GetDumpName("s", ".vtk", step));
  }
  /*
  if (sem("dump")) {
    m.Dump(&as->GetField(), "u");
    m.Dump(&as->GetColorSum(), "cl");
  }
  */
  if (sem("checkloop")) {
    if (as->GetTime() >= tmax) {
      sem.LoopBreak();
    }
    ++step;
  }
  sem.LoopEnd();
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 2
set int by 2
set int bz 1

set int bsx 32
set int bsy 32
set int bsz 1

set int dim 2

set int px 1
set int py 1
set int pz 1

set string init_vf circlels
set vect circle_c 0.5 0.5 0.5
set double circle_r 0.2

set double bcc_clear0 0
set double bcc_clear1 1
set double inletcl 0
set double bcc_fill -1
set string bc_xm symm
set string bc_xp symm
set string bc_ym symm
set string bc_yp symm

set double hypre_symm_tol 1e-6
set double hypre_vort_tol 1e-6
set double hypre_gen_tol 1e-6
set int hypre_periodic_x 1
set int hypre_periodic_y 1
set int hypre_periodic_z 1
)EOF";

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
