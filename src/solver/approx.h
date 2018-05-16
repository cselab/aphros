#include <exception>
#include <memory>
#include <array>

namespace solver {

template <class Idx, class Expr>
class Approx {
 public:
  virtual ~Approx() {}
  virtual Expr GetExpr(Idx) const = 0;
};

template <class M, class Expr>
class InterpolationInnerFaceCentral : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  const M& m;
  explicit InterpolationInnerFaceCentral(const M& m) : m(m) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    Scal a = 0.5;
    e.InsertTerm(1. - a, cm);
    e.InsertTerm(a, cp);
    return e;
  }
};

template <class M, class Expr>
class InterpolationInnerFaceSecondUpwindDeferred 
    : public Approx<IdxFace, Expr> {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const FieldFace<Scal>& probe_;
  const FieldCell<Scal>& fc_prev_;
  const FieldCell<Vect>& fc_prev_grad_;
  const Scal threshold_;
  const M& m;
 public:
  InterpolationInnerFaceSecondUpwindDeferred(
      const M& m,
      const FieldFace<Scal>& probe,
      const FieldCell<Scal>& fc_prev,
      const FieldCell<Vect>& fc_prev_grad,
      Scal threshold = 1e-10)
      : m(m)
      , probe_(probe)
      , fc_prev_(fc_prev)
      , fc_prev_grad_(fc_prev_grad)
      , threshold_(threshold)
  {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    const M& m = this->m;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    // f = fmm + fm + fp
    //const std::array<Scal, 3> a = {0., 1., 0.}; // FOU
    //const std::array<Scal, 3> a = {0., 0.5, 0.5}; // CD
    //const std::array<Scal, 3> a = {-0.5, 1.5, 0.}; // SOU
    const std::array<Scal, 3> a = {-1./8., 6./8., 3./8.}; // QUICK
    if (probe_[f] > threshold_) {
      e.InsertTerm(a[1], cm);
      e.InsertTerm(a[2]+a[0], cp);
      e.SetConstant(4. * a[0]*fc_prev_grad_[cm].dot(m.GetVectToCell(f, 0)));
      //e.InsertTerm(1., cm);
      //e.InsertTerm(0., cp);
      //e.SetConstant(-fc_prev_grad_[cm].dot(m.GetVectToCell(f, 0)));
    } else if (probe_[f] < -threshold_) {
      e.InsertTerm(a[2]+a[0], cm);
      e.InsertTerm(a[1], cp);
      e.SetConstant(4. * a[0]*fc_prev_grad_[cp].dot(m.GetVectToCell(f, 1)));
      //e.InsertTerm(0., cm);
      //e.InsertTerm(1., cp);
      //e.SetConstant(-fc_prev_grad_[cp].dot(m.GetVectToCell(f, 1)));
    } else {
      // TODO: Make a CDS proper for non-uniform grids
      e.InsertTerm(0.5, cm);
      e.InsertTerm(0.5, cp);
    }
    return e;
  }
};

template <class M, class Expr>
class InterpolationBoundaryFaceNearestCell 
    : public Approx<IdxFace, Expr> {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const MapFace<std::shared_ptr<ConditionFace>>& mfc_;
  const M& m;
 public:
  InterpolationBoundaryFaceNearestCell(
      const M& m, const MapFace<std::shared_ptr<ConditionFace>>& mfc)
      : m(m), mfc_(mfc) {}
  Expr GetExpr(IdxFace f) const override {
    const M& m = this->m;
    Expr e;
    if (auto cb = mfc_.find(f)) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      e.InsertTerm(0, cm);
      e.InsertTerm(0, cp);
      if (auto cd = dynamic_cast<ConditionFaceValue<Scal>*>(cb->get())) {
        e.SetConstant(cd->GetValue());
      } else if (auto cd =
          dynamic_cast<ConditionFaceDerivative<Scal>*>(cb->get())) {
        size_t id = cd->GetNci();
        IdxCell cc = m.GetNeighbourCell(f, id);
        Scal g = (id == 0 ? 1. : -1.);
        Scal a = m.GetVectToCell(f, id).norm() * g;
        e.SetConstant(a * cd->GetDerivative());
        e.InsertTerm(1., cc);
      } else {
        throw std::runtime_error("Unknown boundary condition type");
      }
    } else {
      throw std::runtime_error("Boundary condition not set");
    }
    return e;
  }
};

template <class M, class Expr>
class DerivativeInnerFacePlain 
    : public Approx<IdxFace, Expr> {
  using Scal = typename M::Scal;
  const M& m;
 public:
  explicit DerivativeInnerFacePlain(const M& m) : m(m) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    const M& m = this->m;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    auto dm = m.GetVectToCell(f, 0);
    auto dp = m.GetVectToCell(f, 1);
    Scal alpha = Scal(1) / (dp - dm).norm();
    e.InsertTerm(-alpha, cm);
    e.InsertTerm(alpha, cp);
    return e;
  }
};

template <class M, class Expr>
class DerivativeBoundaryFacePlain : public Approx<IdxFace, Expr> {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const MapFace<std::shared_ptr<ConditionFace>>& mfc_;
  const M& m;
 public:
  explicit DerivativeBoundaryFacePlain(
      const M& m, const MapFace<std::shared_ptr<ConditionFace>>& mfc)
      : m(m), mfc_(mfc) {}
  Expr GetExpr(IdxFace f) const override {
    const M& m = this->m;
    Expr e;
    if (auto cb = mfc_.find(f)) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      e.InsertTerm(0, cm);
      e.InsertTerm(0, cp);
      if (auto cd = dynamic_cast<
          ConditionFaceDerivative<Scal>*>(cb->get())) {
        e.SetConstant(cd->GetDerivative());
      } else if (auto cd = dynamic_cast<
          ConditionFaceValue<Scal>*>(cb->get())) {
        size_t id = cd->GetNci();
        IdxCell c = m.GetNeighbourCell(f, id);
        Scal g = (id == 0 ? 1. : -1.);
        Scal a = 1. / m.GetVectToCell(f, id).norm() * g;
        e.SetConstant(a * cd->GetValue());
        e.InsertTerm(-a, c);
      } else {
        throw std::runtime_error("Unknown boundary condition type");
      }
    } else {
      throw std::runtime_error("Boundary condition not set");
    }
    return e;
  }
};

} // namespace solver

