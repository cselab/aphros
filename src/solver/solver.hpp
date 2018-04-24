#pragma once

#include <exception>
#include <memory>

#include "geom/mesh.hpp"
#include "linear/linear.hpp"

template <class Scal>
bool IsNan(Scal a) {
  return !(a * Scal(0) == Scal(0));
}

template <class T, class Idx>
bool IsNan(const GField<T, Idx>& field) {
  for (auto idx : field.GetRange()) {
    if (IsNan(field[idx])) {
      return true;
    }
  }
  return false;
}

template <class T, class Idx>
GMap<T*, Idx> GetPointers(
    const GMap<std::shared_ptr<T>, Idx>& m_ptr_shared) {
  GMap<T*, Idx> m_ptr;
  for (auto it = m_ptr_shared.cbegin(); it != m_ptr_shared.cend(); ++it) {
    m_ptr[it->GetIdx()] = it->GetValue().get();
  }
  return m_ptr;
}

namespace solver {

using IdxCell = IdxCell;
using IdxFace = IdxFace;

class ConditionFace {
 public:
  ConditionFace(size_t nci) : nci_(nci) {}
  virtual ~ConditionFace() {}
  // neighbour cell id
  virtual size_t GetNci() const {
    return nci_;
  }
 private:
  size_t nci_;
};

class ConditionFaceExtrapolation : public ConditionFace {
 public:
  ConditionFaceExtrapolation(size_t nci) : ConditionFace(nci) {}
};

template <class Value>
class ConditionFaceValue : public ConditionFace {
 public:
  ConditionFaceValue(size_t nci) : ConditionFace(nci) {}
  virtual Value GetValue() const = 0;
};

template <class Vect>
class ConditionFaceValueExtractComponent :
    public ConditionFaceValue<typename Vect::value_type> {
  using Scal = typename Vect::value_type;
  using P = ConditionFaceValue<Scal>;
  ConditionFaceValue<Vect>* cond_;
  size_t comp_;
 public:
  ConditionFaceValueExtractComponent(ConditionFaceValue<Vect>* cond, 
                                     size_t comp)
      : P(cond->GetNci())
      , cond_(cond)
      , comp_(comp)
  {}
  Scal GetValue() const override {
    return cond_->GetValue()[comp_];
  }
};

template <class Value>
class ConditionFaceValueFixed : public ConditionFaceValue<Value> {
  Value value_;
 public:
  ConditionFaceValueFixed(const Value& value, size_t nci)
      : ConditionFaceValue<Value>(nci)
      , value_(value)
  {}
  Value GetValue() const override {
    return value_;
  }
  void Set(const Value& v) {
    value_ = v;
  }
};

template <class Value>
class ConditionFaceDerivative : public ConditionFace {
 public:
  ConditionFaceDerivative(size_t nci) : ConditionFace(nci) {}
  virtual Value GetDerivative() const = 0;
};

template <class Value>
class ConditionFaceDerivativeFixed : public ConditionFaceDerivative<Value> {
  Value derivative_;
 public:
  explicit ConditionFaceDerivativeFixed(const Value& derivative, size_t nci)
      : ConditionFaceDerivative<Value>(nci)
      , derivative_(derivative)
  {}
  virtual Value GetDerivative() const override {
    return derivative_;
  }
  void Set(const Value& v) {
    derivative_ = v;
  }
};

class ConditionCell {
 public:
  virtual ~ConditionCell() {}
};

template <class Value>
class ConditionCellValue : public ConditionCell {
 public:
  virtual Value GetValue() const = 0;
};

template <class Value>
class ConditionCellValueFixed : public ConditionCellValue<Value> {
  Value value_;
 public:
  explicit ConditionCellValueFixed(const Value& value)
      : value_(value)
  {}
  Value GetValue() const override {
    return value_;
  }
};

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
class InterpolationInnerFaceFirstUpwind :
    public Approximation<Mesh, IdxFace, Expr> {
  using Scal = typename Mesh::Scal;
  const FieldFace<Scal>& probe_;
  const Scal threshold_;
 public:
  InterpolationInnerFaceFirstUpwind(const Mesh& mesh,
                               const FieldFace<Scal>& probe,
                               Scal threshold = 1e-8)
      : Approximation<Mesh, IdxFace, Expr>(mesh)
      , probe_(probe)
      , threshold_(threshold)
  {}
  Expr GetExpression(IdxFace idxface) const override {
    Expr expr;
    const Mesh& mesh = this->mesh;
    IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
    IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
    if (probe_[idxface] > threshold_) {
      expr.InsertTerm(1., cm);
      expr.InsertTerm(0., cp);
    } else if (probe_[idxface] < -threshold_) {
      expr.InsertTerm(0., cm);
      expr.InsertTerm(1., cp);
    } else {
      expr.InsertTerm(0.5, cm);
      expr.InsertTerm(0.5, cp);
    }
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
    if (probe_[idxface] > threshold_) {
      expr.InsertTerm(1., cm);
      expr.InsertTerm(0., cp);
      expr.SetConstant(-fc_prev_grad_[cm].dot(mesh.GetVectToCell(idxface, 0)));
    } else if (probe_[idxface] < -threshold_) {
      expr.InsertTerm(0., cm);
      expr.InsertTerm(1., cp);
      expr.SetConstant(-fc_prev_grad_[cp].dot(mesh.GetVectToCell(idxface, 1)));
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

template <class T, class Mesh>
FieldFace<T> Interpolate(const FieldNode<T>& fn_u,
                               const Mesh& mesh) {
  FieldFace<T> res(mesh);

  for (auto idxface : mesh.SuFaces()) {
    T sum = static_cast<T>(0);
    for (size_t i = 0; i < mesh.GetNumNeighbourNodes(idxface); ++i) {
      sum += fn_u[mesh.GetNeighbourNode(idxface, i)];
    }
    res[idxface] = sum / mesh.GetNumNeighbourNodes(idxface);
  }

  return res;
}

template <class T, class Mesh>
T GetInterpolatedInner(const FieldCell<T>& fc_u,
                       IdxFace idxface,
                       const Mesh& mesh) {
  using Scal = typename Mesh::Scal;
  IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
  IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
  //Vect xf = mesh.GetCenter(idxface);
  //Vect xm = mesh.GetCenter(cm);
  //Vect xp = mesh.GetCenter(cp);
  //Scal alpha = (xf - xm).dot(xp - xm) / (xp - xm).sqrnorm();
  Scal alpha = 0.5;
  return fc_u[cm] * (1. - alpha) + fc_u[cp] * alpha;
}

template <class Scal>
Scal GetGeometricAverage(Scal a, Scal b) {
  return std::sqrt(a * b);
}

template <class Scal, size_t dim>
GVect<Scal, dim> GetGeometricAverage(
    const GVect<Scal, dim>& a, 
    const GVect<Scal, dim>& b, 
    Scal th = 1e-8) {
  auto m = (a + b) * 0.5;
  if (m.norm() < th) {
    return m;
  }
  auto d = m / m.norm();
  return d * std::sqrt(a.norm() * b.norm());
}

template <class T, class Mesh>
FieldFace<T> Interpolate(
    const FieldCell<T>& fc_u,
    const MapFace<std::shared_ptr<ConditionFace>>& mf_cond_u,
    const Mesh& mesh, bool geometric = false) {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = IdxCell;
  using IdxFace = IdxFace;

  FieldFace<T> res(mesh, T(0)); // Valid value essential for extrapolation

  if (geometric) {
    for (auto idxface : mesh.SuFaces()) {
			IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
			IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
			res[idxface] = GetGeometricAverage(fc_u[cm], fc_u[cp]);
    }
  } else {
    for (auto idxface : mesh.SuFaces()) {
			IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
			IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
			//Vect xf = mesh.GetCenter(idxface);
			//Vect xm = mesh.GetCenter(cm);
			//Vect xp = mesh.GetCenter(cp);
			//Scal alpha = (xf - xm).dot(xp - xm) / (xp - xm).sqrnorm();
			Scal alpha = 0.5;
			res[idxface] = fc_u[cm] * (1. - alpha) + fc_u[cp] * alpha;
    }
  }

  for (auto it = mf_cond_u.cbegin(); it != mf_cond_u.cend(); ++it) {
    IdxFace idxface = it->GetIdx();
    ConditionFace* cond = it->GetValue().get();
    if (auto cond_value = dynamic_cast<ConditionFaceValue<T>*>(cond)) {
      res[idxface] = cond_value->GetValue();
    } else if (auto cond_derivative =
        dynamic_cast<ConditionFaceDerivative<T>*>(cond)) {
      size_t id = cond->GetNci();
      IdxCell cc = mesh.GetNeighbourCell(idxface, id);
      Scal factor = (id == 0 ? 1. : -1.);
      Scal alpha = mesh.GetVectToCell(idxface, id).norm() * factor;
      res[idxface] = fc_u[cc] + cond_derivative->GetDerivative() * alpha;
    } else if (dynamic_cast<ConditionFaceExtrapolation*>(cond)) {
      size_t id = cond->GetNci();
      IdxCell idxcell = mesh.GetNeighbourCell(idxface, id);
      Scal factor = (id == 0 ? 1. : -1.);
      Vect normal = mesh.GetNormal(idxface) * factor;
      Scal dist = mesh.GetVectToCell(idxface, id).norm();
      T nom = fc_u[idxcell] / dist;
      Scal den = 1. / dist;
      Scal volume = mesh.GetVolume(idxcell);
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        IdxFace nface = mesh.GetNeighbourFace(idxcell, i);
        if (nface == idxface) {
          den -= mesh.GetOutwardSurface(idxcell, i).dot(normal) / volume;
        } else {
          nom += res[nface] *
              mesh.GetOutwardSurface(idxcell, i).dot(normal) / volume;
        }
      }
      res[idxface] = nom / den;
    } else {
      throw std::runtime_error("Unknown boundary condition type");
    }
  }
  return res;
}


template <class Mesh>
FieldCell<typename Mesh::Scal>
CalcLaplacian(const FieldCell<typename Mesh::Scal>& fc_u,
              const Mesh& mesh) {
  using Scal = typename Mesh::Scal;
  FieldCell<Scal> res(mesh);

  for (auto idxcell : mesh.Cells()) {
    Scal sum = 0.;
    for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
      auto idxface = mesh.GetNeighbourFace(idxcell, i);
			auto cm = mesh.GetNeighbourCell(idxface, 0);
			auto cp = mesh.GetNeighbourCell(idxface, 1);
			auto dm = mesh.GetVectToCell(idxface, 0);
			auto dp = mesh.GetVectToCell(idxface, 1);
			Scal derivative = (fc_u[cp] - fc_u[cm]) / (dp - dm).norm();
			sum += derivative * mesh.GetArea(idxface) *
					mesh.GetOutwardFactor(idxcell, i);
    }
    res[idxcell] = sum / mesh.GetVolume(idxcell);
  }

  return res;
}

template <class T, class Mesh>
FieldFace<T> InterpolateFirstUpwind(
    const FieldCell<T>& fc_u,
    const MapFace<std::shared_ptr<ConditionFace>>& mf_cond_u,
    const FieldFace<typename Mesh::Scal>& probe,
    const Mesh& mesh, typename Mesh::Scal threshold = 1e-8) {
  using Scal = typename Mesh::Scal;
  using IdxCell = IdxCell;
  using IdxFace = IdxFace;

  FieldFace<T> res(mesh);

  for (IdxFace idxface : mesh.Faces()) {
		 IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
		 IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
		 if (probe[idxface] > threshold) {
			 res[idxface] = fc_u[cm];
		 } else if (probe[idxface] < -threshold) {
			 res[idxface] = fc_u[cp];
		 } else {
			 // TODO: Make a CDS proper for non-uniform grids
			 res[idxface] = 0.5 * (fc_u[cm] + fc_u[cp]);
		 }
  }

  for (auto it = mf_cond_u.cbegin(); it != mf_cond_u.cend(); ++it) {
    IdxFace idxface = it->GetIdx();
    ConditionFace* cond = it->GetValue().get();
    if (auto cond_value = dynamic_cast<ConditionFaceValue<T>*>(cond)) {
      res[idxface] = cond_value->GetValue();
    } else if (auto cond_derivative =
        dynamic_cast<ConditionFaceDerivative<T>*>(cond)) {
      size_t id = cond->GetNci();
      IdxCell cc = mesh.GetNeighbourCell(idxface, id);
      Scal factor = (id == 0 ? 1. : -1.);
      Scal alpha = mesh.GetVectToCell(idxface, id).norm() * factor;
      res[idxface] = fc_u[cc] + cond_derivative->GetDerivative() * alpha;
    } else {
      throw std::runtime_error("Unknown boundary condition type");
    }
  }
  return res;
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

template <class Mesh>
FieldFace<typename Mesh::Scal>
InterpolateSuperbee(
    const FieldCell<typename Mesh::Scal>& fc_u,
    const FieldCell<typename Mesh::Vect>& fc_u_grad,
    const MapFace<std::shared_ptr<ConditionFace>>& mf_cond_u,
    const FieldFace<typename Mesh::Scal>& probe,
    const Mesh& mesh, typename Mesh::Scal threshold = 1e-8) {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = IdxCell;
  using IdxFace = IdxFace;

  FieldFace<Scal> res(mesh);

  for (IdxFace idxface : mesh.SuFaces()) {
		 IdxCell P = mesh.GetNeighbourCell(idxface, 0);
		 IdxCell E = mesh.GetNeighbourCell(idxface, 1);
		 Vect rp = mesh.GetVectToCell(idxface, 0); 
		 Vect re = mesh.GetVectToCell(idxface, 1); 
		 const auto& u = fc_u;
		 const auto& g = fc_u_grad;
		 if (probe[idxface] > threshold) {
			 res[idxface] =
					 u[P] + 0.5 * Superbee(u[E] - u[P],
																 -Scal(4.) * g[P].dot(rp) - (u[E] - u[P]));
		 } else if (probe[idxface] < -threshold) {
			 res[idxface] =
					 u[E] - 0.5 * Superbee(u[E] - u[P],
																 Scal(4.) * g[E].dot(re) - (u[E] - u[P]));
		 } else {
			 // TODO: Make a CDS proper for non-uniform grids
			 res[idxface] = 0.5 * (u[P] + u[E]);
		 }
  }

  // TODO: Move interpolation on boundaries to a function
  for (auto it = mf_cond_u.cbegin(); it != mf_cond_u.cend(); ++it) {
    IdxFace idxface = it->GetIdx();
    ConditionFace* cond = it->GetValue().get();
    if (auto cond_value = dynamic_cast<ConditionFaceValue<Scal>*>(cond)) {
      res[idxface] = cond_value->GetValue();
    } else if (auto cond_derivative =
        dynamic_cast<ConditionFaceDerivative<Scal>*>(cond)) {
      size_t id = cond->GetNci();
      IdxCell cc = mesh.GetNeighbourCell(idxface, id);
      Scal factor = (id == 0 ? 1. : -1.);
      Scal alpha = mesh.GetVectToCell(idxface, id).norm() * factor;
      res[idxface] = fc_u[cc] + cond_derivative->GetDerivative() * alpha;
    } else {
      throw std::runtime_error("Unknown boundary condition type");
    }
  }
  return res;
}

template <class T, class Mesh>
FieldCell<T> Average(const FieldFace<T>& ff_u, const Mesh& mesh) {
  using Scal = typename Mesh::Scal;
  FieldCell<T> res(mesh);
  for (IdxCell idxcell : mesh.AllCells()) {
    T sum(0);
    for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
      IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
      sum += ff_u[idxface];
    }
    res[idxcell] = sum / static_cast<Scal>(mesh.GetNumNeighbourFaces(idxcell));
  }
  return res;
}

template <class T, class Mesh>
void Smoothen(
    FieldCell<T>& fc,
    const MapFace<std::shared_ptr<ConditionFace>>& mf_cond,
    Mesh& m, size_t rep) {
  auto sem = m.GetSem("smoothen");
  for (size_t i = 0; i < rep; ++i) {
    if (sem()) {
      fc = Average(Interpolate(fc, mf_cond, m), m);
      m.Comm(&fc);
    }
  }
}

template <class Mesh>
FieldCell<typename Mesh::Vect> Gradient(
    const FieldFace<typename Mesh::Scal>& ff_u,
    const Mesh& mesh) {
  using Vect = typename Mesh::Vect;
  FieldCell<Vect> res(mesh, Vect::kZero);
  for (auto idxcell : mesh.SuCells()) {
    Vect sum = Vect::kZero;
    for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
      IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
      sum += mesh.GetOutwardSurface(idxcell, i) * ff_u[idxface];
    }
    res[idxcell] = sum / mesh.GetVolume(idxcell);
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

template <class Field, class Mesh, class Scal = typename Mesh::Scal>
Scal CalcDiff(const Field& first, const Field& second, const Mesh& mesh) {
  Scal res = 0.;
  using Idx = typename Field::Idx;
  for (Idx idx : GRange<Idx>(mesh)) {
    res = std::max(res, std::abs(first[idx] - second[idx]));
  }
  return res;
}

template <class Idx, class Mesh, class Scal = typename Mesh::Scal>
Scal CalcDiff(const GField<typename Mesh::Vect, Idx>& first,
                const GField<typename Mesh::Vect, Idx>& second,
                const Mesh& mesh) {
  Scal res = 0.;
  for (Idx idx : mesh.template Get<Idx>()) {
    res = std::max(res, first[idx].dist(second[idx]));
  }
  return res;
}

// derivative(arg_target) = sum_i (coeff_i * func(arg_i))
template <class Scal>
std::vector<Scal> GetDerivativeApproxCoeffs(
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
std::vector<Scal> GetDerivativeApproxCoeffs(
    Scal arg_target, const std::vector<Scal>& args, size_t skip_initial) {
  auto size = args.size();
  size_t cut_size = size - skip_initial;
  std::vector<Scal> cut_args(cut_size);
  for (size_t i = 0; i < cut_size; ++i) {
    cut_args[i] = args[skip_initial + i];
  }
  auto cut_coeffs = GetDerivativeApproxCoeffs(arg_target, cut_args);
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

