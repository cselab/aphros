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

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

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
std::map<Scal, Scal> CalcVolume(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*> fcu,
    const Multi<const FieldCell<Scal>*> fccl, M& m) {
  auto sem = m.GetSem();
  const Scal kClNone = -1;
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
    const Multi<const FieldCell<Scal>*> fck,
    const std::map<Scal, Scal>& mvol0, M& m) {
  auto sem = m.GetSem();
  struct {
    std::map<Scal, Scal> mvol;
  } * ctx(sem);
  auto& mvol = ctx->mvol;

  const Scal kClNone = -1;
  const Scal kvol = 1000.;

  if (sem.Nested("volume")) {
    mvol = CalcVolume(layers, fcu, fccl, m);
  }
  if (sem("calc")) {
    ffv.Reinit(m, 0);
    for (auto f : m.Faces()) {
      IdxCell cm = m.GetCell(f, 0);
      IdxCell cp = m.GetCell(f, 1);
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
        Vect nm(0);
        Vect np(0);
        Scal km = 0;
        Scal kp = 0;
        for (auto i : layers) {
          if ((*fccl[i])[cm] == cl) {
            um = (*fcu[i])[cm];
            nm = (*fcn[i])[cm];
            km = (*fck[i])[cm];
          }
          if ((*fccl[i])[cp] == cl) {
            up = (*fcu[i])[cp];
            np = (*fcn[i])[cp];
            kp = (*fck[i])[cp];
          }
        }
        const Scal k = (std::abs(um - 0.5) < std::abs(up - 0.5) ? km : kp);
        const Vect n = (std::abs(um - 0.5) < std::abs(up - 0.5) ? nm : np);
        const Scal vold =
            mvol0.count(cl)
                ? kvol * ((mvol0.at(cl) - mvol.at(cl)) / mvol0.at(cl))
                : 0;
        const Scal v = m.GetSurface(f).dot(n * (-k + vold));
        if (!IsNan(v)) {
          ffv[f] = v;
        }
      }
    }
  }
}

void Run(M& m, Vars&) {
  auto sem = m.GetSem();

  struct {
    std::unique_ptr<Vofm<M>> as; // advection solver
    FieldCell<Scal> fcs; // volume source
    FieldFace<Scal> ffv; // volume flux
    Multi<FieldCell<Scal>> fck; // curvature
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
  const Scal tmax = 0.1;
  const Scal cfl = 0.5;

  // if (sem.Nested()) {
  //  InitVf(fcu, var, m);
  //}
  if (sem("init")) {
    fcs.Reinit(m, 0);
    ffv.Reinit(m, 0);

    const std::string qpath = "ref/voronoi/cl.dat";
    FieldCell<Scal> qfccl;
    MIdx qsize;
    ReadPlain(qpath, qfccl, qsize);
    GIndex<IdxCell, M::dim> qbc(qsize);
    std::cout << "qsize=" << qbc.GetSize() << std::endl;
    FieldCell<Scal> fccl(m, 0); // initial color
    auto& bc = m.GetIndexCells();
    const MIdx size = m.GetGlobalSize();
    for (auto c : m.Cells()) {
      const MIdx w = bc.GetMIdx(c);
      const MIdx qw = w * qsize / size;
      const IdxCell qc = qbc.GetIdx(qw);
      fccl[c] = qfccl[qc];
    }

    FieldCell<Scal> fcu(m, 1); // initial volume fraction

    typename Vofm<M>::Par p;
    p.sharpen = true;
    p.sharpen_cfl = 0.1;
    p.recolor_unionfind = true;
    p.recolor_reduce = true;
    p.recolor_grid = true;
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
        ffv, layers, as->GetFieldM(), as->GetColor(), as->GetNormal(), fck,
        mvol0, m);
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
    as->SetTimeStep(dt);
    if (m.IsRoot()) {
      std::cout << "dt=" << dt << std::endl;
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
  if (sem("dump")) {
    m.Dump(&as->GetField(), "u");
    m.Dump(&as->GetColorSum(), "cl");
  }
  if (sem("checkloop0")) {
  }
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
set int bx 8
set int by 8
set int bz 1

set int bsx 8
set int bsy 8
set int bsz 1

set int dim 2

set int px 1
set int py 1
set int pz 1

set string init_vf circlels
set vect circle_c 0.5 0.5 0.5
set double circle_r 0.2
)EOF";

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
