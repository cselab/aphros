#include <exception>
#include <memory>
#include <array>

namespace solver {

template <class M, class Idx, class Expr>
class Approximation {
  using Scal = typename M::Scal;
 protected:
  const M& m;
 public:
  explicit Approximation(const M& m)
      : m(m)
  {}
  virtual Expr GetExpression(Idx idx) const = 0;
};

template <class M, class Expr>
class InterpolationInnerFaceCentral :
    public Approximation<M, IdxFace, Expr> {
  using Scal = typename M::Scal;
 public:
  explicit InterpolationInnerFaceCentral(const M& m)
      : Approximation<M, IdxFace, Expr>(m)
  {}
  Expr GetExpression(IdxFace f) const override {
    Expr expr;
    const M& m = this->m;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    //auto xf = m.GetCenter(f);
    //auto xm = m.GetCenter(cm);
    //auto xp = m.GetCenter(cp);
    //Scal alpha = (xf - xm).dot(xp - xm) / (xp - xm).sqrnorm();
    Scal alpha = 0.5;
    expr.InsertTerm(1. - alpha, cm);
    expr.InsertTerm(alpha, cp);
    return expr;
  }
};

template <class M, class Expr>
class InterpolationInnerFaceSecondUpwindDeferred :
    public Approximation<M, IdxFace, Expr> {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const FieldFace<Scal>& probe_;
  const FieldCell<Scal>& fc_prev_;
  const FieldCell<Vect>& fc_prev_grad_;
  const Scal threshold_;
 public:
  InterpolationInnerFaceSecondUpwindDeferred(
      const M& m,
      const FieldFace<Scal>& probe,
      const FieldCell<Scal>& fc_prev,
      const FieldCell<Vect>& fc_prev_grad,
      Scal threshold = 1e-10)
      : Approximation<M, IdxFace, Expr>(m)
      , probe_(probe)
      , fc_prev_(fc_prev)
      , fc_prev_grad_(fc_prev_grad)
      , threshold_(threshold)
  {}
  Expr GetExpression(IdxFace f) const override {
    Expr expr;
    const M& m = this->m;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    // f = fmm + fm + fp
    //const std::array<Scal, 3> a = {0., 1., 0.}; // FOU
    //const std::array<Scal, 3> a = {0., 0.5, 0.5}; // CD
    //const std::array<Scal, 3> a = {-0.5, 1.5, 0.}; // SOU
    const std::array<Scal, 3> a = {-1./8., 6./8., 3./8.}; // QUICK
    if (probe_[f] > threshold_) {
      expr.InsertTerm(a[1], cm);
      expr.InsertTerm(a[2]+a[0], cp);
      expr.SetConstant(4. * a[0]*fc_prev_grad_[cm].dot(m.GetVectToCell(f, 0)));
      //expr.InsertTerm(1., cm);
      //expr.InsertTerm(0., cp);
      //expr.SetConstant(-fc_prev_grad_[cm].dot(m.GetVectToCell(f, 0)));
    } else if (probe_[f] < -threshold_) {
      expr.InsertTerm(a[2]+a[0], cm);
      expr.InsertTerm(a[1], cp);
      expr.SetConstant(4. * a[0]*fc_prev_grad_[cp].dot(m.GetVectToCell(f, 1)));
      //expr.InsertTerm(0., cm);
      //expr.InsertTerm(1., cp);
      //expr.SetConstant(-fc_prev_grad_[cp].dot(m.GetVectToCell(f, 1)));
    } else {
      // TODO: Make a CDS proper for non-uniform grids
      expr.InsertTerm(0.5, cm);
      expr.InsertTerm(0.5, cp);
    }
    return expr;
  }
};

template <class M, class Expr>
class InterpolationBoundaryFaceNearestCell:
    public Approximation<M, IdxFace, Expr> {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const MapFace<std::shared_ptr<ConditionFace>>& mf_cond_;
 public:
  InterpolationBoundaryFaceNearestCell(
      const M& m,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_cond)
      : Approximation<M, IdxFace, Expr>(m)
      , mf_cond_(mf_cond)
  {}
  Expr GetExpression(IdxFace f) const override {
    const M& m = this->m;
    Expr expr;
    if (auto cond_generic = mf_cond_.find(f)) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      expr.InsertTerm(0, cm);
      expr.InsertTerm(0, cp);
      if (auto cond_value =
          dynamic_cast<ConditionFaceValue<Scal>*>(cond_generic->get())) {
        expr.SetConstant(cond_value->GetValue());
      } else if (auto cond_derivative =
          dynamic_cast<ConditionFaceDerivative<Scal>*>(cond_generic->get())) {
        size_t id = cond_derivative->GetNci();
        IdxCell cc = m.GetNeighbourCell(f, id);
        Scal factor = (id == 0 ? 1. : -1.);
        Scal alpha = m.GetVectToCell(f, id).norm() * factor;
        expr.SetConstant(alpha * cond_derivative->GetDerivative());
        expr.InsertTerm(1., cc);
      } else {
        throw std::runtime_error("Unknown boundary condition type");
      }
    } else {
      throw std::runtime_error("Boundary condition not set");
    }
    return expr;
  }
};

template <class M, class Expr>
class DerivativeInnerFacePlain:
    public Approximation<M, IdxFace, Expr> {
  using Scal = typename M::Scal;
 public:
  explicit DerivativeInnerFacePlain(const M& m)
      : Approximation<M, IdxFace, Expr>(m)
  {}
  Expr GetExpression(IdxFace f) const override {
    Expr expr;
    const M& m = this->m;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    auto dm = m.GetVectToCell(f, 0);
    auto dp = m.GetVectToCell(f, 1);
    Scal alpha = Scal(1) / (dp - dm).norm();
    expr.InsertTerm(-alpha, cm);
    expr.InsertTerm(alpha, cp);
    return expr;
  }
};

template <class M, class Expr>
class DerivativeBoundaryFacePlain:
    public Approximation<M, IdxFace, Expr> {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const MapFace<std::shared_ptr<ConditionFace>>& mf_cond_;
 public:
  explicit DerivativeBoundaryFacePlain(
      const M& m,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_cond)
      : Approximation<M, IdxFace, Expr>(m)
      , mf_cond_(mf_cond)
  {}
  Expr GetExpression(IdxFace f) const override {
    const M& m = this->m;
    Expr expr;
    if (auto cond_generic = mf_cond_.find(f)) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      expr.InsertTerm(0, cm);
      expr.InsertTerm(0, cp);
      if (auto cond_derivative =
          dynamic_cast<ConditionFaceDerivative<Scal>*>(cond_generic->get())) {
        expr.SetConstant(cond_derivative->GetDerivative());
      } else if (auto cond_value =
          dynamic_cast<ConditionFaceValue<Scal>*>(cond_generic->get())) {
        size_t id = cond_value->GetNci();
        IdxCell cc = m.GetNeighbourCell(f, id);
        Scal factor = (id == 0 ? 1. : -1.);
        Scal alpha = 1. / m.GetVectToCell(f, id).norm() * factor;
        expr.SetConstant(alpha * cond_value->GetValue());
        expr.InsertTerm(-alpha, cc);
      } else {
        throw std::runtime_error("Unknown boundary condition type");
      }
    } else {
      throw std::runtime_error("Boundary condition not set");
    }
    return expr;
  }
};

} // namespace solver

