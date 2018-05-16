#pragma once

#include <exception>
#include <memory>

#include "geom/mesh.h"
#include "linear/linear.h"
#include "cond.h"
#include "approx.h"

template <class Scal>
bool IsNan(Scal a) {
  return !(a * Scal(0) == Scal(0));
}

template <class T, class Idx>
bool IsNan(const GField<T, Idx>& u) {
  for (auto i : u.GetRange()) {
    if (IsNan(u[i])) {
      return true;
    }
  }
  return false;
}

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
  using Vect = typename M::Vect;

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
  using Vect = typename M::Vect;

  for (auto f : m.SuFaces()) {
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    Scal a = 0.5;
    ff[f] = fc[cm] * (1. - a) + fc[cp] * a;
  }
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
    if (auto cd = dynamic_cast<CondFaceVal<T>*>(cb)) {
      ff[f] = cd->GetValue();
    } else if (auto cd = dynamic_cast<CondFaceGrad<T>*>(cb)) {
      size_t id = cb->GetNci();
      IdxCell c = m.GetNeighbourCell(f, id);
      Scal w = (id == 0 ? 1. : -1.);
      Scal a = m.GetVectToCell(f, id).norm() * w;
      ff[f] = fc[c] + cd->GetGrad() * a;
    } else if (dynamic_cast<CondFaceExtrap*>(cb)) {
      // TODO test
      size_t id = cb->GetNci();
      IdxCell c = m.GetNeighbourCell(f, id);
      Scal w = (id == 0 ? 1. : -1.);
      Vect n = m.GetNormal(f) * w;
      Scal h = m.GetVectToCell(f, id).norm();
      T a = fc[c] / h;
      Scal b = 1. / h;
      Scal vol = m.GetVolume(c);
      for (auto q : m.Nci(c)) {
        IdxFace fq = m.GetNeighbourFace(c, q);
        if (fq == f) {
          b -= m.GetOutwardSurface(c, q).dot(n) / vol;
        } else {
          a += ff[fq] * m.GetOutwardSurface(c, q).dot(n) / vol;
        }
      }
      ff[f] = a / b;
    } else {
      // TODO add name to CondFace etc
      throw std::runtime_error("Unknown boundary condition type");
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
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  FieldFace<T> ff(m, T(0)); // Valid 0 needed for CondFaceExtrap

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
// fc: field cell [a]
// fcg: gradient of field [a]
// mfc: face cond
// ffw: flow direction [s]
// Output:
// field face [s]
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
			 ff[f] = u[cp] + 0.5 * Superbee(du, -4. * g[cm].dot(rm) - du);
		 } else if (ffw[f] < -th) {
			 ff[f] = u[cm] - 0.5 * Superbee(du, 4. * g[cp].dot(rp) - du);
		 } else {
			 ff[f] = 0.5 * (u[cm] + u[cp]);
		 }
  }

  InterpolateB(fc, mfc, ff, m);
  return ff;
}

template <class T, class M>
FieldCell<T> Average(const FieldFace<T>& ff_u, const M& m) {
  using Scal = typename M::Scal;
  FieldCell<T> res(m);
  for (IdxCell c : m.AllCells()) {
    T sum(0);
    for (size_t i = 0; i < m.GetNumNeighbourFaces(c); ++i) {
      IdxFace f = m.GetNeighbourFace(c, i);
      sum += ff_u[f];
    }
    res[c] = sum / static_cast<Scal>(m.GetNumNeighbourFaces(c));
  }
  return res;
}

template <class T, class M>
void Smoothen(
    FieldCell<T>& fc,
    const MapFace<std::shared_ptr<CondFace>>& mf_cond,
    M& m, size_t rep) {
  auto sem = m.GetSem("smoothen");
  for (size_t i = 0; i < rep; ++i) {
    if (sem()) {
      fc = Average(Interpolate(fc, mf_cond, m), m);
      m.Comm(&fc);
    }
  }
}

template <class M>
FieldCell<typename M::Vect> Gradient(
    const FieldFace<typename M::Scal>& ff_u,
    const M& m) {
  using Vect = typename M::Vect;
  FieldCell<Vect> res(m, Vect::kZero);
  for (auto c : m.SuCells()) {
    Vect sum = Vect::kZero;
    for (size_t i = 0; i < m.GetNumNeighbourFaces(c); ++i) {
      IdxFace f = m.GetNeighbourFace(c, i);
      sum += m.GetOutwardSurface(c, i) * ff_u[f];
    }
    res[c] = sum / m.GetVolume(c);
  }
  return res;
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
  size_t iter_;

 protected:
  void IncIter() {
    ++iter_;
  }
  void ClearIter() {
    iter_ = 0;
  }

 public:
  UnsteadyIterativeSolver(double t, double dt)
      : UnsteadySolver(t, dt)
      , iter_(0)
  {}
  virtual void MakeIteration() = 0;
  virtual double GetError() const {
    return 0.;
  }
  virtual size_t GetIter() const {
    return iter_;
  }
  virtual void StartStep() override {
    ClearIter();
  }
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

template <class T>
struct LayersData {
  T time_curr, time_prev, iter_curr, iter_prev;
  T& Get(Layers layer) {
    switch (layer) {
      case Layers::time_curr: {
        return time_curr;
      }
      case Layers::time_prev: {
        return time_prev;
      }
      case Layers::iter_curr: {
        return iter_curr;
      }
      case Layers::iter_prev: {
        return iter_prev;
      }
      default: {
        throw(std::runtime_error("LayersData::Get(): Unknown layer"));
      }
    }
  }
  const T& Get(Layers layer) const {
    return const_cast<const T&>(const_cast<LayersData*>(this)->Get(layer));
  }
};

// Convention: Use Get/Set for fast procedures and Calc for those requiring computation:
// GetValue(field, idx) vs GetNorm(field)

template <class Field, class M, class Scal = typename M::Scal>
Scal CalcDiff(const Field& first, const Field& second, const M& m) {
  Scal res = 0.;
  using Idx = typename Field::Idx;
  for (Idx idx : GRange<Idx>(m)) {
    res = std::max(res, std::abs(first[idx] - second[idx]));
  }
  return res;
}

template <class Idx, class M, class Scal = typename M::Scal>
Scal CalcDiff(const GField<typename M::Vect, Idx>& first,
                const GField<typename M::Vect, Idx>& second,
                const M& m) {
  Scal res = 0.;
  for (Idx idx : m.template Get<Idx>()) {
    res = std::max(res, first[idx].dist(second[idx]));
  }
  return res;
}

// derivative(arg_target) = sum_i (coeff_i * func(arg_i))
template <class Scal>
std::vector<Scal> GetGradCoeffs(
    Scal arg_target, const std::vector<Scal>& args) {

  auto size = args.size();
  std::vector<Scal> res(size);
  for (size_t i = 0; i < size; ++i) {
    Scal denom = 1.;
    Scal numer = 0.;
    for (size_t j = 0; j < size; ++j) {
      if (j != i) {
        denom *= args[i] - args[j];
        Scal term = 1.;
        for (size_t k = 0; k < size; ++k) {
          if (k != i && k != j) {
            term *= arg_target - args[k];
          }
        }
        numer += term;
      }
    }
    res[i] = numer / denom;
  }
  return res;
}


template <class Scal>
std::vector<Scal> GetGradCoeffs(
    Scal arg_target, const std::vector<Scal>& args, size_t skip_initial) {
  auto size = args.size();
  size_t cut_size = size - skip_initial;
  std::vector<Scal> cut_args(cut_size);
  for (size_t i = 0; i < cut_size; ++i) {
    cut_args[i] = args[skip_initial + i];
  }
  auto cut_coeffs = GetGradCoeffs(arg_target, cut_args);
  std::vector<Scal> res(size);
  for (size_t i = 0; i < skip_initial; ++i) {
    res[i] = 0.;
  }
  for (size_t i = 0; i < cut_size; ++i) {
    res[skip_initial + i] = cut_coeffs[i];
  }
  return res;
}

} // namespace solver

