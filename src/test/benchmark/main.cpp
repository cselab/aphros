// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>

#include "geom/mesh.h"
#include "geom/rangemulti.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/solver.h"
#include "util/format.h"
#include "util/sysinfo.h"
#include "util/timer.h"

#define EXPOSE(x)                 \
  do {                            \
    volatile auto EXPOSE_vol = x; \
    (void)EXPOSE_vol;             \
  } while (0);

const int dim = 3;
using MIdx = generic::MIdx<dim>;
using IdxCell = IdxCell;
using IdxFace = IdxFace;
using Dir = GDir<dim>;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;
using M = MeshStructured<Scal, dim>;

// Echo Execute
#define EE(...)                                   \
  ;                                               \
  std::cerr << "\n" << #__VA_ARGS__ << std::endl; \
  __VA_ARGS__;

M GetMesh(MIdx size) {
  const Rect<Vect> dom(Vect(0), Vect(1));
  const MIdx begin(0);
  const int halos = 2;
  return InitUniformMesh<M>(dom, begin, size, halos, true, true, size, 0);
}

// Covered range of cells
enum class Cover { all, su, in };

class TimerMesh : public ExecutionTimer {
 public:
  TimerMesh(const std::string& name, M& m_, Cover cover = Cover::in)
      : ExecutionTimer(name, 0.05, 10), m(m_), cover_(cover) {}
  Cover GetCover() const {
    return cover_;
  }

 protected:
  M& m;
  Cover cover_;
};

class LoopAllCellsPlain : public TimerMesh {
 public:
  LoopAllCellsPlain(M& m_) : TimerMesh("loop-allcells-plain", m_, Cover::all) {}
  void F() override {
    size_t a = 0;
    for (size_t i = 0; i < m.GetAllBlockCells().size(); ++i) {
      a += i;
    }
    EXPOSE(a);
  }
};

class LoopInCellsPlain : public TimerMesh {
 public:
  LoopInCellsPlain(M& m_) : TimerMesh("loop-incells-plain", m_) {}
  void F() override {
    size_t a = 0;
    for (size_t i = 0; i < m.GetInBlockCells().size(); ++i) {
      a += i;
    }
    EXPOSE(a);
  }
};

class LoopAllCells : public TimerMesh {
 public:
  LoopAllCells(M& m_) : TimerMesh("loop-allcells", m_, Cover::all) {}
  void F() override {
    size_t a = 0;
    for (auto i : m.AllCells()) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

class LoopInCells : public TimerMesh {
 public:
  LoopInCells(M& m_) : TimerMesh("loop-incells", m_) {}
  void F() override {
    size_t a = 0;
    for (auto i : m.Cells()) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

class LoopInCellsNew : public TimerMesh {
 public:
  LoopInCellsNew(M& m_) : TimerMesh("loop-incells-new", m_) {}
  void F() override {
    size_t a = 0;
    const auto indexc = m.GetIndexCells();
    const auto blockci = m.GetInBlockCells();
    using RangeMulti = typename generic::RangeMulti<IdxCell, dim, 1>;
    using MIdxR = typename RangeMulti::MIdx;
    const MIdxR size(blockci.GetSize());
    const MIdxR lead(indexc.GetSize());
    const RangeMulti range(indexc.GetIdx(blockci.GetBegin()), size, lead);
    for (auto i : range) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

class LoopAllFaces : public TimerMesh {
 public:
  LoopAllFaces(M& m_) : TimerMesh("loop-allfaces", m_, Cover::all) {}
  void F() override {
    size_t a = 0;
    for (auto i : m.AllFaces()) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

class LoopInFaces : public TimerMesh {
 public:
  LoopInFaces(M& m_) : TimerMesh("loop-infaces", m_) {}
  void F() override {
    size_t a = 0;
    for (auto i : m.Faces()) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

// Loop access to field
class LoopFldPlain : public TimerMesh {
 public:
  LoopFldPlain(M& m_)
      : TimerMesh("loop-fld-plain", m_, Cover::all)
      , v(GRange<IdxCell>(m).size()) {}
  void F() override {
    Scal a = 0;
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] += a;
      a = v[i];
    }
    EXPOSE(a);
  }

 private:
  std::vector<Scal> v;
};

class LoopFldAllCells : public TimerMesh {
 public:
  LoopFldAllCells(M& m_)
      : TimerMesh("loop-fld-allcells", m_, Cover::all), v(m) {}
  void F() override {
    Scal a = 0;
    for (auto i : m.AllCells()) {
      v[i] += a;
      a = v[i];
    }
    EXPOSE(a);
  }

 private:
  FieldCell<Scal> v;
};

class LoopMIdxAllCells : public TimerMesh {
 public:
  LoopMIdxAllCells(M& m_) : TimerMesh("loop-midx-allcells", m_, Cover::all) {}
  void F() override {
    size_t a = 0;
    for (auto c : m.AllCells()) {
      auto w = m.GetIndexCells().GetMIdx(c);
      a += w[0];
    }
    EXPOSE(a);
  }
};

class LoopMIdxAllFaces : public TimerMesh {
 public:
  LoopMIdxAllFaces(M& m_) : TimerMesh("loop-midx-allfaces", m_, Cover::all) {}
  void F() override {
    size_t a = 0;
    for (auto f : m.AllFaces()) {
      auto wd = m.GetIndexFaces().GetMIdxDir(f);
      a += wd.first[0];
      a += wd.second.GetLetter();
    }
    EXPOSE(a);
  }
};

class CellVolume : public TimerMesh {
 public:
  CellVolume(M& m_) : TimerMesh("m-cell-volume", m_) {}
  void F() override {
    size_t a = 1;
    for (auto c : m.Cells()) {
      a += m.GetVolume(c);
    }
    EXPOSE(a);
  }
};

class CellCenter : public TimerMesh {
 public:
  CellCenter(M& m_) : TimerMesh("m-cell-center", m_) {}
  void F() override {
    int a = 0;
    for (auto c : m.Cells()) {
      a += m.GetCenter(c)[0];
    }
    EXPOSE(a);
  }
};

class CellNCell : public TimerMesh {
 public:
  CellNCell(M& m_) : TimerMesh("m-cell-n-cell", m_) {}
  void F() override {
    size_t a = 1;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetCell(c, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class CellNFace : public TimerMesh {
 public:
  CellNFace(M& m_) : TimerMesh("m-cell-n-face", m_) {}
  void F() override {
    size_t a = 1;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetFace(c, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class CellOutward : public TimerMesh {
 public:
  CellOutward(M& m_) : TimerMesh("m-cell-outward", m_) {}
  void F() override {
    size_t a = 0;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetOutwardFactor(c, q);
      }
    }
    EXPOSE(a);
  }
};

class CellNNode : public TimerMesh {
 public:
  CellNNode(M& m_) : TimerMesh("m-cell-n-node", m_) {}
  void F() override {
    size_t a = 1;
    for (auto c : m.Cells()) {
      for (size_t q = 0; q < m.GetNumNodes(c); ++q) {
        a += m.GetNode(c, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class FaceCenter : public TimerMesh {
 public:
  FaceCenter(M& m_) : TimerMesh("m-face-center", m_) {}
  void F() override {
    int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetCenter(f)[0];
    }
    EXPOSE(a);
  }
};

class FaceSurf : public TimerMesh {
 public:
  FaceSurf(M& m_) : TimerMesh("m-face-surf", m_) {}
  void F() override {
    int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetSurface(f)[0];
    }
    EXPOSE(a);
  }
};

class FaceArea : public TimerMesh {
 public:
  FaceArea(M& m_) : TimerMesh("m-face-area", m_) {}
  void F() override {
    int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetSurface(f)[0];
    }
    EXPOSE(a);
  }
};

class FaceNCell : public TimerMesh {
 public:
  FaceNCell(M& m_) : TimerMesh("m-face-n-cell", m_) {}
  void F() override {
    size_t a = 1;
    for (auto f : m.Faces()) {
      for (size_t q = 0; q < m.GetNumCells(f); ++q) {
        a += m.GetCell(f, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class FaceNNode : public TimerMesh {
 public:
  FaceNNode(M& m_) : TimerMesh("m-face-n-node", m_) {}
  void F() override {
    size_t a = 1;
    for (auto f : m.Faces()) {
      for (size_t q = 0; q < m.GetNumNodes(f); ++q) {
        a += m.GetNode(f, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class Interp : public TimerMesh {
 public:
  Interp(M& m_) : TimerMesh("interp", m_), fc(m), ff(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
    auto& bf = m.GetIndexFaces();
    for (auto f : m.Faces()) {
      if (bf.GetMIdx(f)[0] == 0 && bf.GetDir(f) == Dir(0)) {
        mfc[f].type = BCondType::neumann;
      }
    }
    assert(mfc.size() > 0);
  }
  void F() override {
    size_t a = 0;
    ff = UEmbed<M>::Interpolate(fc, mfc, m);
    a = ff[IdxFace(a)];
    EXPOSE(a);
  }

 private:
  FieldCell<Scal> fc;
  MapEmbed<BCond<Scal>> mfc;
  FieldFace<Scal> ff;
};

class Grad : public TimerMesh {
 public:
  Grad(M& m_) : TimerMesh("grad", m_), fc(m), ff(m) {
    for (auto i : m.AllFaces()) {
      ff[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    static size_t a = 0;
    fc = UEmbed<M>::Gradient(ff, m);
    a += fc[IdxCell(0)][0];
    EXPOSE(a);
  }

 private:
  FieldCell<Vect> fc;
  FieldFace<Scal> ff;
};

class GradTemplate : public TimerMesh {
 public:
  GradTemplate(M& m_) : TimerMesh("grad-template", m_), fc(m), ff(m) {
    for (auto i : m.AllFaces()) {
      ff[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    static size_t a = 0;
    m.ForEachCell([&](auto c) {
      Vect sum(0);
      for (auto q : m.Nci(c)) {
        auto f = c.face(q);
        sum[f.direction] += ff[f] * c.outward_factor(q) / c.h[f.direction];
      }
      fc[c] = sum;
    });
    a += fc[IdxCell(0)][0];
    EXPOSE(a);
  }

 private:
  FieldCell<Vect> fc;
  FieldFace<Scal> ff;
};

class ExplVisc : public TimerMesh {
 public:
  ExplVisc(M& m_)
      : TimerMesh("explvisc", m_, Cover::su), fcv(m), fcf(m), ffmu(m) {
    for (auto i : m.AllCells()) {
      auto a = i.GetRaw();
      fcv[i] = Vect(std::sin(a), std::sin(a + 1), std::sin(a + 2));
    }
    for (auto i : m.AllFaces()) {
      auto a = i.GetRaw();
      ffmu[i] = std::sin(a);
    }
  }
  void F() override {
    static size_t a = 0;
    for (size_t n = 0; n < dim; ++n) {
      FieldCell<Scal> fc = GetComponent(fcv, n);
      auto ff = UEmbed<M>::Interpolate(fc, {}, m);
      auto gc = UEmbed<M>::Gradient(ff, m);
      auto gf = UEmbed<M>::Interpolate(gc, {}, m);
      for (auto c : m.SuCells()) {
        Vect sum(0);
        for (auto q : m.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          sum += gf[f] * (ffmu[f] * m.GetOutwardSurface(c, q)[n]);
        }
        fcf[c] += sum / m.GetVolume(c);
      }
    }
    a += fcf[IdxCell(0)][0];
    EXPOSE(a);
  }

 private:
  FieldCell<Vect> fcv;
  FieldCell<Vect> fcf;
  FieldFace<Scal> ffmu;
};

class ExplViscNotation : public TimerMesh {
 public:
  ExplViscNotation(M& m_)
      : TimerMesh("explvisc-notation", m_, Cover::su), fcv(m), fcf(m), ffmu(m) {
    for (auto i : m.AllCells()) {
      auto a = i.GetRaw();
      fcv[i] = Vect(std::sin(a), std::sin(a + 1), std::sin(a + 2));
    }
    for (auto i : m.AllFaces()) {
      auto a = i.GetRaw();
      ffmu[i] = std::sin(a);
    }
  }
  void F() override {
    static size_t a = 0;
    for (size_t n = 0; n < dim; ++n) {
      FieldCell<Scal> fc = GetComponent(fcv, n);
      auto ff = UEmbed<M>::Interpolate(fc, {}, m);
      auto gc = UEmbed<M>::Gradient(ff, m);
      auto gf = UEmbed<M>::Interpolate(gc, {}, m);
      for (auto c : m.SuCellsM()) {
        Vect sum(0);
        for (auto q : m.Nci(c)) {
          const auto f = c.face(q);
          sum += gf[f] * (ffmu[f] * c.outward_surface(q)[n]);
        }
        fcf[c] += sum / c.volume;
      }
    }
    a += fcf[IdxCell(0)][0];
    EXPOSE(a);
  }

 private:
  FieldCell<Vect> fcv;
  FieldCell<Vect> fcf;
  FieldFace<Scal> ffmu;
};

template <class Field>
class GradientOperatorCell {
 public:
  GradientOperatorCell(const Field& ff_, const M& m_) : ff(ff_), m(m_) {}
  Vect operator[](typename M::Cell c) const {
    Vect res(0);
    for (auto q : m.Nci(c)) {
      auto f = c.face(q);
      res[f.direction] += ff[f] * c.outward_factor(q) / c.h[f.direction];
    }
    return res;
  }

  const Field& ff;
  const M& m;
};

template <class Field>
class GradientOperator {
 public:
  GradientOperator(
      const Field& fc_, const MapEmbed<BCond<Scal>>& mbc_, const M& m_)
      : fc(fc_), mebc(mbc_), m(m_) {}
  auto operator[](typename M::Face f) const {
    return (fc[f.cell(1)] - fc[f.cell(0)]) / f.h[0];
  }

  const Field& fc;
  const MapEmbed<BCond<Scal>>& mebc;
  const M& m;
};

template <class Field>
class InterpolateOperator {
 public:
  InterpolateOperator(
      const Field& fc_, const MapEmbed<BCond<Scal>>& mbc_, const M& m_)
      : fc(fc_), mebc(mbc_), m(m_) {}
  auto operator[](typename M::Face f) const {
    return (fc[f.cell(1)] + fc[f.cell(0)]) * 0.5;
  }

  const Field& fc;
  const MapEmbed<BCond<Scal>>& mebc;
  const M& m;
};

template <class Field>
class ComponentOperator {
 public:
  ComponentOperator(const Field& fc_, size_t d_) : fc(fc_), d(d_) {}
  auto operator[](typename M::Cell c) const {
    return fc[c][d];
  }

  const Field& fc;
  const size_t d;
};

namespace func {

template <class Field>
auto Gradient(const Field& fc, const MapEmbed<BCond<Scal>>& bc, const M& m) {
  return GradientOperator<Field>(fc, bc, m);
}

template <class Field>
auto Interpolate(const Field& fc, const MapEmbed<BCond<Scal>>& bc, const M& m) {
  return InterpolateOperator<Field>(fc, bc, m);
}

template <class Field>
auto Gradient(const Field& ff, const M& m) {
  return GradientOperatorCell<Field>(ff, m);
}

template <class Field>
auto Component(const Field& fc, size_t d) {
  return ComponentOperator<Field>(fc, d);
}

} // namespace func

class ExplViscTemplate : public TimerMesh {
 public:
  ExplViscTemplate(M& m_)
      : TimerMesh("explvisc-template", m_, Cover::su), fcv(m), fcf(m), ffmu(m) {
    for (auto i : m.AllCells()) {
      auto a = i.GetRaw();
      fcv[i] = Vect(std::sin(a), std::sin(a + 1), std::sin(a + 2));
    }
    for (auto i : m.AllFaces()) {
      auto a = i.GetRaw();
      ffmu[i] = std::sin(a);
    }
  }
  void F() override {
    static size_t a = 0;
    for (size_t n = 0; n < dim; ++n) {
      auto fc = GetComponent(fcv, n);
      auto ff = func::Interpolate(fc, {}, m);
      auto gc = func::Gradient(ff, m);
      auto gf = func::Interpolate(gc, {}, m);
      m.ForEachCell([&](auto c) { //
        Vect sum(0);
        for (auto q : m.Nci(c)) {
          const auto f = c.face(q);
          sum += gf[f] * (ffmu[f] * c.outward_factor(q) *
                          (n == f.direction ? f.area : 0));
        }
        fcf[c] += sum / c.volume;
      });
    }
    a += fcf[IdxCell(0)][0];
    EXPOSE(a);
  }

 private:
  FieldCell<Vect> fcv;
  FieldCell<Vect> fcf;
  FieldFace<Scal> ffmu;
};

class Stencil : public TimerMesh {
 public:
  Stencil(M& m_) : TimerMesh("stencil", m_) {}
  void F() override {
    size_t a = 0;
    for (auto c : m.Cells()) {
      for (auto cn : m.Stencil(c)) {
        a += size_t(cn);
      }
    }
    EXPOSE(a);
  }
};

class StencilField : public TimerMesh {
 public:
  StencilField(M& m_) : TimerMesh("stencil_field", m_), fc(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    size_t a = 0;
    for (auto c : m.Cells()) {
      for (auto cn : m.Stencil(c)) {
        a += fc[cn];
      }
    }
    EXPOSE(a);
  }

 private:
  FieldCell<Scal> fc;
};

// test: test index
// m: mesh
// Output:
// time: total per one call [sec]
// iters: number of calls
// mem: memory usage in bytes
// cover: covered range of cells
// name: test name
// Returns 1 if test with index test found
bool RunTest(
    const size_t test, M& m, /*out*/ double& time, size_t& iters, size_t& mem,
    Cover& cover, std::string& name) {
  size_t k = 0;
  TimerMesh* p = nullptr;

  auto create = [&](auto* kernel) {
    if (k++ == test) {
      using T = typename std::remove_pointer<decltype(kernel)>::type;
      p = new T(m);
    }
  };

  create((LoopAllCellsPlain*)0);
  create((LoopInCellsPlain*)0);
  create((LoopAllCells*)0);
  create((LoopInCells*)0);
  create((LoopInCellsNew*)0);
  create((LoopAllFaces*)0);
  create((LoopInFaces*)0);
  create((LoopFldPlain*)0);
  create((LoopFldAllCells*)0);
  create((LoopMIdxAllCells*)0);
  create((LoopMIdxAllFaces*)0);
  create((Interp*)0);
  create((Grad*)0);
  create((GradTemplate*)0);
  create((ExplVisc*)0);
  create((ExplViscNotation*)0);
  create((ExplViscTemplate*)0);

  create((CellVolume*)0);
  create((CellCenter*)0);
  create((FaceCenter*)0);
  create((FaceSurf*)0);
  create((FaceArea*)0);
  create((CellNCell*)0);
  create((CellNFace*)0);
  create((CellOutward*)0);
  create((CellNNode*)0);
  create((FaceNCell*)0);
  create((FaceNNode*)0);

  create((Stencil*)0);
  create((StencilField*)0);

  if (!p) {
    return false;
  }

  auto e = p->Run();
  time = e.min_call_time;
  iters = e.iters;
  mem = sysinfo::GetMem();
  cover = p->GetCover();
  name = p->GetName();
  delete p;

  return true;
}

int main() {
  const std::vector<MIdx> meshsizes = {MIdx(8), MIdx(16), MIdx(32)};

  const std::string fmt = "{:22}{:13.2f}{:8}{:8}{:10.1f}{:17}\n";

  std::stringstream header;
  using std::setw;
  header << util::Format(
      fmt, //
      "name", "t/cell[ns]", "cover", "iters", "mem[MB]", "mem/allcells[B]");

  for (auto meshsize : meshsizes) {
    size_t mem0 = sysinfo::GetMem();
    auto m = GetMesh(meshsize);
    const size_t allcells = m.GetAllBlockCells().size();
    const size_t sucells = m.GetSuBlockCells().size();
    const size_t incells = m.GetInBlockCells().size();
    std::cout << "Mesh " << meshsize << " allcells=" << allcells
              << " incells=" << incells << std::endl;
    std::cout << header.str();

    int test = 0;
    double time;
    size_t iters;
    size_t mem;
    Cover cover;
    std::string name;
    while (RunTest(test++, m, time, iters, mem, cover, name)) {
      size_t dmem = mem - mem0;
      const auto covcells =
          (cover == Cover::all  ? allcells
           : cover == Cover::in ? incells
                                : sucells);
      const auto covname =
          (cover == Cover::all  ? "all"
           : cover == Cover::in ? "in"
                                : "su");
      std::cout << util::Format(
          fmt, //
          name, time * 1e9 / covcells, covname, iters, (dmem / double(1 << 20)),
          (dmem / allcells));
    }
    std::cout << std::endl;
  }
}
