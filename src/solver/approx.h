#include <exception>
#include <memory>
#include <array>

namespace solver {

template <class Mesh, class Idx, class Expr>
class Approximation {
  using Scal = typename Mesh::Scal;
 protected:
  const Mesh& mesh;
 public:
  explicit Approximation(const Mesh& mesh)
      : mesh(mesh)
  {}
  virtual Expr GetExpression(Idx idx) const = 0;
};

template <class Mesh, class Expr>
class InterpolationInnerFaceCentral :
    public Approximation<Mesh, IdxFace, Expr> {
  using Scal = typename Mesh::Scal;
 public:
  explicit InterpolationInnerFaceCentral(const Mesh& mesh)
      : Approximation<Mesh, IdxFace, Expr>(mesh)
  {}
  Expr GetExpression(IdxFace idxface) const override {
    Expr expr;
    const Mesh& mesh = this->mesh;
    IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
    IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
    //auto xf = mesh.GetCenter(idxface);
    //auto xm = mesh.GetCenter(cm);
    //auto xp = mesh.GetCenter(cp);
    //Scal alpha = (xf - xm).dot(xp - xm) / (xp - xm).sqrnorm();
    Scal alpha = 0.5;
    expr.InsertTerm(1. - alpha, cm);
    expr.InsertTerm(alpha, cp);
    return expr;
  }
};

template <class Mesh, class Expr>
class InterpolationInnerFaceSecondUpwindDeferred :
    public Approximation<Mesh, IdxFace, Expr> {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  const FieldFace<Scal>& probe_;
  const FieldCell<Scal>& fc_prev_;
  const FieldCell<Vect>& fc_prev_grad_;
  const Scal threshold_;
 public:
  InterpolationInnerFaceSecondUpwindDeferred(
      const Mesh& mesh,
      const FieldFace<Scal>& probe,
      const FieldCell<Scal>& fc_prev,
      const FieldCell<Vect>& fc_prev_grad,
      Scal threshold = 1e-10)
      : Approximation<Mesh, IdxFace, Expr>(mesh)
      , probe_(probe)
      , fc_prev_(fc_prev)
      , fc_prev_grad_(fc_prev_grad)
      , threshold_(threshold)
  {}
  Expr GetExpression(IdxFace idxface) const override {
    Expr expr;
    const Mesh& mesh = this->mesh;
    IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
    IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
    // f = fmm + fm + fp
    //const std::array<Scal, 3> a = {0., 1., 0.}; // FOU
    const std::array<Scal, 3> a = {0., 0.5, 0.5}; // CD
    //const std::array<Scal, 3> a = {-0.5, 1.5, 0.}; // SOU
    //const std::array<Scal, 3> a = {-1./8., 6./8., 3./8.}; // QUICK
    if (probe_[idxface] > threshold_) {
      expr.InsertTerm(a[1], cm);
      expr.InsertTerm(a[2]+a[0], cp);
      expr.SetConstant(4. * a[0]*fc_prev_grad_[cm].dot(mesh.GetVectToCell(idxface, 0)));
      //expr.InsertTerm(1., cm);
      //expr.InsertTerm(0., cp);
      //expr.SetConstant(-fc_prev_grad_[cm].dot(mesh.GetVectToCell(idxface, 0)));
    } else if (probe_[idxface] < -threshold_) {
      expr.InsertTerm(a[2]+a[0], cm);
      expr.InsertTerm(a[1], cp);
      expr.SetConstant(4. * a[0]*fc_prev_grad_[cp].dot(mesh.GetVectToCell(idxface, 1)));
      //expr.InsertTerm(0., cm);
      //expr.InsertTerm(1., cp);
      //expr.SetConstant(-fc_prev_grad_[cp].dot(mesh.GetVectToCell(idxface, 1)));
    } else {
      // TODO: Make a CDS proper for non-uniform grids
      expr.InsertTerm(0.5, cm);
      expr.InsertTerm(0.5, cp);
    }
    return expr;
  }
};

template <class Mesh, class Expr>
class InterpolationBoundaryFaceNearestCell:
    public Approximation<Mesh, IdxFace, Expr> {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  const MapFace<std::shared_ptr<ConditionFace>>& mf_cond_;
 public:
  InterpolationBoundaryFaceNearestCell(
      const Mesh& mesh,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_cond)
      : Approximation<Mesh, IdxFace, Expr>(mesh)
      , mf_cond_(mf_cond)
  {}
  Expr GetExpression(IdxFace idxface) const override {
    const Mesh& mesh = this->mesh;
    Expr expr;
    if (auto cond_generic = mf_cond_.find(idxface)) {
      IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
      IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
      expr.InsertTerm(0, cm);
      expr.InsertTerm(0, cp);
      if (auto cond_value =
          dynamic_cast<ConditionFaceValue<Scal>*>(cond_generic->get())) {
        expr.SetConstant(cond_value->GetValue());
      } else if (auto cond_derivative =
          dynamic_cast<ConditionFaceDerivative<Scal>*>(cond_generic->get())) {
        size_t id = cond_derivative->GetNci();
        IdxCell cc = mesh.GetNeighbourCell(idxface, id);
        Scal factor = (id == 0 ? 1. : -1.);
        Scal alpha = mesh.GetVectToCell(idxface, id).norm() * factor;
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

template <class Mesh, class Expr>
class DerivativeInnerFacePlain:
    public Approximation<Mesh, IdxFace, Expr> {
  using Scal = typename Mesh::Scal;
 public:
  explicit DerivativeInnerFacePlain(const Mesh& mesh)
      : Approximation<Mesh, IdxFace, Expr>(mesh)
  {}
  Expr GetExpression(IdxFace idxface) const override {
    Expr expr;
    const Mesh& mesh = this->mesh;
    IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
    IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
    auto dm = mesh.GetVectToCell(idxface, 0);
    auto dp = mesh.GetVectToCell(idxface, 1);
    Scal alpha = Scal(1) / (dp - dm).norm();
    expr.InsertTerm(-alpha, cm);
    expr.InsertTerm(alpha, cp);
    return expr;
  }
};

template <class Mesh, class Expr>
class DerivativeBoundaryFacePlain:
    public Approximation<Mesh, IdxFace, Expr> {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  const MapFace<std::shared_ptr<ConditionFace>>& mf_cond_;
 public:
  explicit DerivativeBoundaryFacePlain(
      const Mesh& mesh,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_cond)
      : Approximation<Mesh, IdxFace, Expr>(mesh)
      , mf_cond_(mf_cond)
  {}
  Expr GetExpression(IdxFace idxface) const override {
    const Mesh& mesh = this->mesh;
    Expr expr;
    if (auto cond_generic = mf_cond_.find(idxface)) {
      IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
      IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
      expr.InsertTerm(0, cm);
      expr.InsertTerm(0, cp);
      if (auto cond_derivative =
          dynamic_cast<ConditionFaceDerivative<Scal>*>(cond_generic->get())) {
        expr.SetConstant(cond_derivative->GetDerivative());
      } else if (auto cond_value =
          dynamic_cast<ConditionFaceValue<Scal>*>(cond_generic->get())) {
        size_t id = cond_value->GetNci();
        IdxCell cc = mesh.GetNeighbourCell(idxface, id);
        Scal factor = (id == 0 ? 1. : -1.);
        Scal alpha = 1. / mesh.GetVectToCell(idxface, id).norm() * factor;
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

