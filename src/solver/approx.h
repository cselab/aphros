#pragma once

#include <exception>
#include <memory>
#include <array>
#include <stdexcept>
#include <string>

namespace solver {

template <class Idx, class Expr>
class Approx {
 public:
  virtual ~Approx() {}
  virtual Expr GetExpr(Idx) const = 0;
};

// Convection scheme
enum class ConvSc { fou, cd, sou, quick };

// Convection scheme name
std::string GetName(ConvSc sc);

// Convection scheme by name.
ConvSc GetConvSc(std::string s);

// sc: scheme:
//   - fou: first order upwind
//   - cd: central differences (mean value)
//   - sou: second order upwind
//   - quick: QUICK
template <class Scal>
std::array<Scal, 3> GetCoeff(ConvSc sc) {
  std::array<Scal, 3> a;
  switch (sc) {
    case ConvSc::fou: 
      a = {0., 1., 0.}; 
      break; 
    case ConvSc::cd: 
      a = {0., 0.5, 0.5}; 
      break; 
    case ConvSc::sou: 
      a = {-0.5, 1.5, 0.}; 
      break;
    case ConvSc::quick: 
      a = {-1./8., 6./8., 3./8.}; 
      break;
    default:
      throw std::runtime_error("GetCoeff: invalid ConvSc");
  }
  return a;
}

// Explicit interpolation to faces near inner cells.
// fc: field cell [s]
// fc: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M>
void InterpolateI(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const FieldFace<T>& ffw,  const M& m, ConvSc sc, typename M::Scal th,
    FieldFace<T>& ff) {
  using Scal = typename M::Scal;

  // f = fmm*a[0] + fm*a[1] + fp*a[2]
  std::array<Scal, 3> a = GetCoeff<Scal>(sc);

  ff.Reinit(m);
  for (auto f : m.Faces()) {
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    if (ffw[f] > th) {
      ff[f] = 4. * a[0] * fcgp[cm].dot(m.GetVectToCell(f, 0)) +
          a[1] * fc[cm] + (a[2] + a[0]) * fc[cp];
    } else if (ffw[f] < -th) {
      ff[f] = 4. * a[0] * fcgp[cp].dot(m.GetVectToCell(f, 1)) +
          a[1] * fc[cp] + (a[2] + a[0]) * fc[cm];
    } else {
      ff[f] = (fc[cm] + fc[cp]) * 0.5;
    }
  }
}

// Interpolates from cells to inner faces.
// T: value type (Scal or Vect)
// fc: field cell [s]
// fcgp: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M>
void Interpolate(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const MapFace<std::shared_ptr<CondFace>>& mfc, 
    const FieldFace<T>& ffw,  const M& m, ConvSc sc, typename M::Scal th,
    FieldFace<T>& ff) {
  InterpolateI(fc, fcgp, ffw, m, sc, th, ff);
  InterpolateB(fc, mfc, ff, m);
}

// Implicit interpolation to inner faces with deferred correction.
// fc: field cell [s]
// fc: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
//   - fou: first order upwind
//   - cd: central differences (mean value)
//   - sou: second order upwind
//   - quick: QUICK
// df: deferred correction factor (1. fully deferred)
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M, class Expr>
void InterpolateI(const FieldCell<T>& fc,
                  const FieldCell<typename M::Vect>& fcgp,
                  const FieldFace<T>& ffw, FieldFace<Expr>& ff, 
                  const M& m, ConvSc sc, typename M::Scal df,
                  typename M::Scal th) {
  using Scal = typename M::Scal;

  // f = fmm*a[0] + fm*a[1] + fp*a[2]
  std::array<Scal, 3> a = GetCoeff<Scal>(sc);

  ff.Reinit(m);
  auto dfm = 1. - df;
  for (auto f : m.Faces()) {
    Expr& e = ff[f];
    e.Clear();
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    if (ffw[f] > th) {
      e.InsertTerm(a[1] * dfm + df, cm);
      e.InsertTerm((a[2] + a[0]) * dfm, cp);
      e.SetConstant(4. * a[0] * fcgp[cm].dot(m.GetVectToCell(f, 0)) +
          (a[1] - 1.) * df * fc[cm] + (a[2] + a[0]) * df * fc[cp]);
    } else if (ffw[f] < -th) {
      e.InsertTerm((a[2] + a[0]) * dfm, cm);
      e.InsertTerm(a[1] * dfm + df, cp);
      e.SetConstant(4. * a[0] * fcgp[cp].dot(m.GetVectToCell(f, 1)) +
          (a[1] - 1.) * df * fc[cp] + (a[2] + a[0]) * df * fc[cm]);
    } else {
      e.InsertTerm(0.5, cm);
      e.InsertTerm(0.5, cp);
    }
  }
}

template <class M, class Expr>
class FaceValB : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  FaceValB(const M& m, const MapFace<std::shared_ptr<CondFace>>& mfc)
      : m(m), mfc_(mfc) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    if (auto cb = mfc_.find(f)) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      e.InsertTerm(0, cm);
      e.InsertTerm(0, cp);
      if (auto cd = dynamic_cast<CondFaceVal<Scal>*>(cb->get())) {
        e.SetConstant(cd->GetValue());
      } else if (auto cd = dynamic_cast<CondFaceGrad<Scal>*>(cb->get())) {
        size_t id = cd->GetNci();
        IdxCell c = m.GetNeighbourCell(f, id);
        Scal g = (id == 0 ? 1. : -1.);
        Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * g;
        e.SetConstant(a * cd->GetGrad());
        e.InsertTerm(1., c);
      } else {
        throw std::runtime_error("FaceValB: unknown cond");
      }
    } else {
      throw std::runtime_error("FaceValB: unset cond");
    }
    return e;
  }

 private:
  const M& m;
  const MapFace<std::shared_ptr<CondFace>>& mfc_;
};

template <class M, class Expr>
class FaceGrad : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  explicit FaceGrad(const M& m) : m(m) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    // XXX: adhoc uniform
    Scal a = m.GetArea(f) / m.GetVolume(cp);
    e.InsertTerm(-a, cm);
    e.InsertTerm(a, cp);
    return e;
  }

 private:
  const M& m;
};

// Explicit gradient on faces near inner cells.
// fc: field [s]
// Output:
// ff: normal gradient [i]
template <class M, class T>
void GradientI(const FieldCell<T>& fc, const M& m, FieldFace<T>& ff) {
  using Scal = typename M::Scal;

  ff.Reinit(m);
  for (auto f : m.Faces()) {
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    Scal a = m.GetArea(f) / m.GetVolume(cp);
    ff[f] = (fc[cp] - fc[cm]) * a;
  }
}

// Explicit gradient on boundary faces.
// fc: field [s]
// mfc: face conditions
// Output:
// ff: normal gradient [i]
template <class M, class T>
void GradientB(const FieldCell<T>& fc,
    const MapFace<std::shared_ptr<CondFace>>& mfc,
    const M& m, FieldFace<T>& ff) {
  using Scal = typename M::Scal;

  for (const auto& it : mfc) {
    IdxFace f = it.GetIdx();
    CondFace* cb = it.GetValue().get(); // cond base
    if (auto cd = dynamic_cast<CondFaceGrad<T>*>(cb)) {
      ff[f] = cd->GetGrad();
    } else if (auto cd = dynamic_cast<CondFaceVal<T>*>(cb)) {
      size_t id = cd->GetNci();
      IdxCell c = m.GetNeighbourCell(f, id);
      Scal g = (id == 0 ? 1. : -1.);
      Scal hr = m.GetArea(f) / m.GetVolume(c);
      Scal a = hr * 2 * g;
      ff[f] = (cd->GetValue() - fc[c]) * a;
    } else {
      throw std::runtime_error("GradientB: unknown cond");
    }
  }
}

// Explicit gradient on inner faces.
// fc: field [s]
// Output:
// ff: normal gradient [i]
template <class M, class T>
void Gradient(
    const FieldCell<T>& fc, const MapFace<std::shared_ptr<CondFace>>& mfc,
    const M& m, FieldFace<T>& ff) {
  GradientI(fc, m, ff);
  GradientB(fc, mfc, m, ff);
}

// Implicit gradient in inner faces.
// Output:
// ff: normal gradient [i]
template <class M, class Expr>
void GradientI(FieldFace<Expr>& ff, const M& m) {
  using Scal = typename M::Scal;

  ff.Reinit(m);
  for (auto f : m.Faces()) {
    Expr& e = ff[f];
    e.Clear();

    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    // XXX: adhoc uniform
    Scal a = m.GetArea(f) / m.GetVolume(cp);
    e.InsertTerm(-a, cm);
    e.InsertTerm(a, cp);
  }
}


template <class M, class Expr>
class FaceGradB : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  explicit FaceGradB(
      const M& m, const MapFace<std::shared_ptr<CondFace>>& mfc)
      : m(m), mfc_(mfc) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    if (auto cb = mfc_.find(f)) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      e.InsertTerm(0, cm);
      e.InsertTerm(0, cp);
      if (auto cd = dynamic_cast<CondFaceGrad<Scal>*>(cb->get())) {
        e.SetConstant(cd->GetGrad());
      } else if (auto cd = dynamic_cast<CondFaceVal<Scal>*>(cb->get())) {
        size_t id = cd->GetNci();
        IdxCell c = m.GetNeighbourCell(f, id);
        Scal g = (id == 0 ? 1. : -1.);
        Scal hr = m.GetArea(f) / m.GetVolume(c);
        Scal a = hr * 2 * g;
        e.SetConstant(a * cd->GetValue());
        e.InsertTerm(-a, c);
      } else {
        throw std::runtime_error("FaceGradB: unknown cond");
      }
    } else {
      throw std::runtime_error("FaceGradB: unset cond");
    }
    return e;
  }

 private:
  const M& m;
  const MapFace<std::shared_ptr<CondFace>>& mfc_;
};

} // namespace solver

