#pragma once

#include <exception>
#include <memory>
#include <cmath>
#include <string>

#include "geom/mesh.h"
#include "cond.h"
#include "approx.h"


namespace solver {

// Interpolates from nodes to faces
template <class T, class M>
FieldFace<T> Interpolate(const FieldNode<T>& fn, const M& m) {
  FieldFace<T> ff(m);
  for (auto f : m.SuFaces()) {
    T s(0);
    for (size_t i = 0; i < m.GetNumNeighbourNodes(f); ++i) {
      s += fn[m.GetNeighbourNode(f, i)];
    }
    ff[f] = s / m.GetNumNeighbourNodes(f);
  }
  return ff;
}

// Interpolation to inner faces.
// fc: field cell [s]
// Output:
// ff: face cell [i]
template <class T, class M>
void InterpolateI(const FieldCell<T>& fc, FieldFace<T>& ff, const M& m) {
  using Scal = typename M::Scal;

  for (auto f : m.Faces()) {
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    Scal a = 0.5;
    ff[f] = fc[cm] * (1. - a) + fc[cp] * a;
  }
}

// Interpolation to support faces.
// fc: field cell [a]
// Output:
// ff: face cell [s]
template <class T, class M>
void InterpolateS(const FieldCell<T>& fc, FieldFace<T>& ff, const M& m) {
  using Scal = typename M::Scal;

  for (auto f : m.SuFaces()) {
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    Scal a = 0.5;
    ff[f] = fc[cm] * (1. - a) + fc[cp] * a;
  }
}

template <class Scal>
class UReflectFace {
 public:
  using Vect = GVect<Scal, 3>;
  // v: value
  // n: normal to face
  static Scal Get(Scal v, const Vect& /*n*/) {
    return v;
  }
  static Vect Get(const Vect& v, const Vect& n) {
    return v - n * n.dot(v);
  }
};

template <class Scal>
class UReflectCell {
 public:
  using Vect = GVect<Scal, 3>;
  // v: value
  // n: normal to face
  static Scal Get(Scal v, const Vect& /*n*/) {
    return v;
  }
  static Vect Get(const Vect& v, const Vect& n) {
    return v - n * (2. * n.dot(v));
  }
};


// Linear extrapolation.
// xt: target
// x0,x1: points
// v0,v1: values
template <class T, class Scal>
T UExtrap(Scal xt, Scal x0, const T& v0, Scal x1, const T& v1) {
  return v0 + (v1 - v0) * ((xt - x0) / (x1 - x0));
}

// Interpolation to faces with defined conditions.
// fc: field cell [i]
// mfc: face cond
// Output:
// ff: values updated on faces defined in mfc
template <class T, class M>
void InterpolateB(
    const FieldCell<T>& fc,
    const MapFace<std::shared_ptr<CondFace>>& mfc, 
    FieldFace<T>& ff, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  for (const auto& it : mfc) {
    IdxFace f = it.GetIdx();
    CondFace* cb = it.GetValue().get(); // cond base
    size_t nci = cb->GetNci();
    if (auto cd = dynamic_cast<CondFaceVal<T>*>(cb)) {
      ff[f] = cd->GetValue();
    } else if (auto cd = dynamic_cast<CondFaceGrad<T>*>(cb)) {
      IdxCell c = m.GetNeighbourCell(f, nci);
      Scal w = (nci == 0 ? 1. : -1.);
      Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * w;
      ff[f] = fc[c] + cd->GetGrad() * a;
    } else if (dynamic_cast<CondFaceExtrap*>(cb)) {
      // TODO test
      IdxCell c = m.GetNeighbourCell(f, nci);
      size_t q = m.GetNci(c, f);
      size_t qo = m.GetOpposite(q);
      IdxFace fo = m.GetNeighbourFace(c, qo);
      Vect n = m.GetNormal(f);
      // cell 
      const T& v0 = fc[c];
      Scal x0 = 0.;
      // opposite face
      const T& v1 = ff[fo];
      Scal x1 = n.dot(m.GetCenter(fo) - m.GetCenter(c));
      // target 
      Scal xt = n.dot(m.GetCenter(f) - m.GetCenter(c));
      
      ff[f] = UExtrap(xt, x0, v0, x1, v1);
    } else if (dynamic_cast<CondFaceReflect*>(cb)) {
      // TODO test
      IdxCell c = m.GetNeighbourCell(f, nci);
      Vect n = m.GetNormal(f);
      auto v = fc[c];
      ff[f] = UReflectFace<Scal>::Get(v, n);
    } else {
      // TODO add name to CondFace etc
      throw std::runtime_error("InterpolateB: unknown cond");
    }
  }
}

// Interpolates from cells to support faces.
// T: value type (Scal or Vect)
// fc: field cell [a]
// mfc: face cond
// Output:
// field face [s]
template <class T, class M>
FieldFace<T> Interpolate(
    const FieldCell<T>& fc,
    const MapFace<std::shared_ptr<CondFace>>& mfc, 
    const M& m) {
  FieldFace<T> ff(m); // Valid 0 needed for CondFaceExtrap

  InterpolateS(fc, ff, m);
  InterpolateB(fc, mfc, ff, m);

  return ff;
}

template <class Scal>
Scal Superbee(Scal p, Scal q) {
  if(p > 0. && q > 0.) {
    return std::max(std::min(2*p, q), std::min(p, 2*q));
  } else if(p < 0. && q < 0.) {
    return -std::max(std::min(-2*p, -q), std::min(-p, -2*q));
  }
  return 0.;
}

// Second order upwind interpolation with TVD Superbee limiter
// fc: fieldcell [a]
// fcg: gradient of field [a]
// mfc: face cond
// ffw: flow direction [s]
// Output:
// fieldface [s]
template <class M>
FieldFace<typename M::Scal> InterpolateSuperbee(
    const FieldCell<typename M::Scal>& fc,
    const FieldCell<typename M::Vect>& fcg,
    const MapFace<std::shared_ptr<CondFace>>& mfc,
    const FieldFace<typename M::Scal>& ffw,
    const M& m, typename M::Scal th = 1e-8) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  FieldFace<Scal> ff(m);

  for (IdxFace f : m.SuFaces()) {
		 IdxCell cm = m.GetNeighbourCell(f, 0);
		 IdxCell cp = m.GetNeighbourCell(f, 1);
		 Vect rm = m.GetVectToCell(f, 0); 
		 Vect rp = m.GetVectToCell(f, 1); 
		 const auto& u = fc;
		 const auto& g = fcg;
     Scal du = u[cp] - u[cm];
		 if (ffw[f] > th) {
			 ff[f] = u[cm] + 0.5 * Superbee(du, -4. * g[cm].dot(rm) - du);
		 } else if (ffw[f] < -th) {
			 ff[f] = u[cp] - 0.5 * Superbee(du, 4. * g[cp].dot(rp) - du);
		 } else {
			 ff[f] = 0.5 * (u[cm] + u[cp]);
		 }
  }

  InterpolateB(fc, mfc, ff, m);
  return ff;
}

// Returns average of fieldface.
// ff: fieldface [a]
// Output:
// fieldcell [a]
template <class T, class M>
FieldCell<T> Average(const FieldFace<T>& ff, const M& m) {
  using Scal = typename M::Scal;
  FieldCell<T> fc(m);
  for (IdxCell c : m.AllCells()) {
    T s(0);
    for (auto q : m.Nci(c)) {
      IdxFace f = m.GetNeighbourFace(c, q);
      s += ff[f];
    }
    fc[c] = s / Scal(m.GetNumNeighbourFaces(c));
  }
  return fc;
}

// Smoothens fieldcell.
// fc: fieldcell [s]
// mfc: condface
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void Smoothen(FieldCell<T>& fc,
              const MapFace<std::shared_ptr<CondFace>>& mfc,
              M& m, size_t rep) {
  auto sem = m.GetSem("smoothen");
  for (size_t i = 0; i < rep; ++i) {
    if (sem()) {
      fc = Average(Interpolate(fc, mfc, m), m);
      m.Comm(&fc);
    }
  }
}

// Returns gradient.
// ff: scalar fieldface [s]
// Output:
// gradient vector [s]
template <class M>
FieldCell<typename M::Vect> Gradient(
    const FieldFace<typename M::Scal>& ff, const M& m) {
  using Vect = typename M::Vect;
  FieldCell<Vect> fc(m);
  for (auto c : m.SuCells()) {
    Vect s(0);
    for (auto q : m.Nci(c)) {
      IdxFace f = m.GetNeighbourFace(c, q);
      s += m.GetOutwardSurface(c, q) * ff[f];
    }
    fc[c] = s / m.GetVolume(c);
  }
  return fc;
}

class UnsteadySolver {
 public:
  UnsteadySolver(double t, double dt) : t_(t), dt_(dt) {}
  virtual ~UnsteadySolver() {}
  virtual void StartStep() {}
  virtual void FinishStep() { IncTime(); }
  virtual double GetTime() const { return t_; }
  virtual void SetTime(double t) { t_  = t; }
  virtual double GetTimeStep() const { return dt_; }
  virtual void SetTimeStep(double dt) { dt_ = dt; }

 protected:
  virtual void IncTime() { t_ += dt_; }

 private:
  double t_;
  double dt_;
};

class UnsteadyIterativeSolver : public UnsteadySolver {
 public:
  UnsteadyIterativeSolver(double t, double dt)
      : UnsteadySolver(t, dt) , i_(0) {}
  virtual void MakeIteration() = 0;
  virtual double GetError() const { return 0.; }
  virtual size_t GetIter() const { return i_; }
  virtual void StartStep() override { ClearIter(); }

 protected:
  void IncIter() { ++i_; }
  void ClearIter() { i_ = 0; }

 private:
  size_t i_;
};

// Between StartStep() and FinishStep():
//   iter_curr -- last iteration of next step
//   iter_prev -- previous iteration of next step
//   time_curr -- current step
//   time_prev -- previous step
// After FinishStep() call:
//   iter_curr, iter_prev -- undefined
//   time_curr -- current step
//   time_prev -- previous step
enum class Layers { time_curr, time_prev, iter_curr, iter_prev };

std::string GetName(Layers);

template <class T>
struct LayersData {
  T time_curr, time_prev, iter_curr, iter_prev;
  const T& Get(Layers l) const {
    switch (l) {
      case Layers::time_curr: { return time_curr; }
      case Layers::time_prev: { return time_prev; }
      case Layers::iter_curr: { return iter_curr; }
      case Layers::iter_prev: { return iter_prev; }
      default: {
        throw std::runtime_error("LayersData::Get(): Unknown layer");
      }
    }
  }
  T& Get(Layers l) {
    return const_cast<T&>(const_cast<const LayersData*>(this)->Get(l));
  }
};

// Convention: Use Get/Set for fast procedures and Calc for those requiring computation:
// GetValue(field, idx) vs GetNorm(field)

template <class Field, class M, class Scal = typename M::Scal>
Scal CalcDiff(const Field& fa, const Field& fb, const M& m) {
  Scal r = 0.;
  using Idx = typename Field::Idx;
  for (Idx i : m.template GetIn<Idx>()) {
    r = std::max<Scal>(r, std::abs(fa[i] - fb[i]));
  }
  return r;
}

template <class Idx, class M, class Scal = typename M::Scal>
Scal CalcDiff(const GField<typename M::Vect, Idx>& fa,
                const GField<typename M::Vect, Idx>& fb,
                const M& m) {
  Scal r = 0.;
  for (Idx i : m.template GetIn<Idx>()) {
    r = std::max<Scal>(r, (fa[i] - fb[i]).norminf());
  }
  return r;
}

// Coefficients for approximation of gradient with polynomial.
// x: target point
// z: stencil points
// Output:
// k: such that grad(x) = sum_i (ki * f(zi))
template <class Scal>
std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z) {
  // TODO: test

  size_t s = z.size();
  std::vector<Scal> k(s);
  for (size_t i = 0; i < s; ++i) {
    Scal a = 0.;
    Scal b = 1.;
    for (size_t j = 0; j < s; ++j) {
      if (j != i) {
        b *= z[i] - z[j];
        Scal t = 1.;
        for (size_t k = 0; k < s; ++k) {
          if (k != i && k != j) {
            t *= x - z[k];
          }
        }
        a += t;
      }
    }
    k[i] = a / b;
  }
  return k;
}


// Returns GetGradCoeffs(x,z[b:]) preceeded by b zeros.
template <class Scal>
std::vector<Scal> GetGradCoeffs(
    Scal x, const std::vector<Scal>& z, size_t b) {
  size_t s = z.size();
  size_t ss = s - b;
  std::vector<Scal> zz(ss);
  for (size_t i = 0; i < ss; ++i) {
    zz[i] = z[b + i];
  }
  std::vector<Scal> kk = GetGradCoeffs(x, zz);
  std::vector<Scal> k(s);
  for (size_t i = 0; i < b; ++i) {
    k[i] = 0.;
  }
  for (size_t i = 0; i < ss; ++i) {
    k[b + i] = kk[i];
  }
  return k;
}

// apply reflection to field on Reflect boundaries or other if force=true
// fill: value for other types that CondFaceReflect
template <class T, class M>
void BcReflect(FieldCell<T>& uc,
               const MapFace<std::shared_ptr<CondFace>>& mfc,
               T fill, bool force, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Dir = typename M::Dir;
  auto& bf = m.GetIndexFaces();
  auto& bc = m.GetIndexCells();
  for (const auto& it : mfc) {
    CondFace* cb = it.GetValue().get();
    size_t nci = cb->GetNci();
    IdxFace f = it.GetIdx();
    Dir df = bf.GetDir(f);
    Vect n = m.GetNormal(f);
    // offset from face towards cell (inner normal to boundary)
    MIdx wo(0);
    wo[size_t(df)] = (nci == 0 ? -1 : 1);
    IdxCell cp = m.GetNeighbourCell(f, nci);
    MIdx wp = bc.GetMIdx(cp);
    MIdx wpp = wp + wo;
    MIdx wm = wp - wo;
    MIdx wmm = wm - wo;
    IdxCell cm = bc.GetIdx(wm);
    IdxCell cmm = bc.GetIdx(wmm);
    IdxCell cpp = bc.GetIdx(wpp);
    // apply
    if (dynamic_cast<CondFaceReflect*>(cb) || force) {
      uc[cm] = UReflectCell<Scal>::Get(uc[cp], n);
      uc[cmm] = UReflectCell<Scal>::Get(uc[cpp], n);
    } else {
      uc[cm] = fill;
      uc[cmm] = fill;
    }
  }
}

} // namespace solver

