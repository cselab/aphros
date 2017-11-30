/*
 * solver.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "mesh.hpp"
#include "linear.hpp"
#include <exception>
#include <memory>

using geom::IntIdx;

template <class Scal>
bool IsNan(Scal a) {
  return !(a * Scal(0) == Scal(0));
}

template <class T, class Idx>
bool IsNan(const geom::FieldGeneric<T, Idx>& field) {
  for (auto idx : field.GetRange()) {
    if (IsNan(field[idx])) {
      return true;
    }
  }
  return false;
}

template <class T, class Idx>
geom::MapGeneric<T*, Idx> GetPointers(
    const geom::MapGeneric<std::shared_ptr<T>, Idx>& m_ptr_shared) {
  geom::MapGeneric<T*, Idx> m_ptr;
  for (auto it = m_ptr_shared.cbegin(); it != m_ptr_shared.cend(); ++it) {
    m_ptr[it->GetIdx()] = it->GetValue().get();
  }
  return m_ptr;
}

namespace solver {

using IdxCell = geom::IdxCell;
using IdxFace = geom::IdxFace;

class ConditionFace {
 public:
  virtual ~ConditionFace() {}
};

class ConditionFaceExtrapolation : public ConditionFace {
 public:
  ConditionFaceExtrapolation() {}
};

template <class ValueType>
class ConditionFaceValue : public ConditionFace {
 public:
  virtual ValueType GetValue() const = 0;
};

template <class Vect>
class ConditionFaceValueExtractComponent :
    public ConditionFaceValue<typename Vect::value_type> {
  using Scal = typename Vect::value_type;
  ConditionFaceValue<Vect>* cond_;
  size_t comp_;
 public:
  ConditionFaceValueExtractComponent(ConditionFaceValue<Vect>* cond, size_t comp)
      : cond_(cond)
      , comp_(comp)
  {}
  Scal GetValue() const override {
    return cond_->GetValue()[comp_];
  }
};

template <class ValueType>
class ConditionFaceValueFixed : public ConditionFaceValue<ValueType> {
  ValueType value_;
 public:
  ConditionFaceValueFixed(const ValueType& value)
      : value_(value)
  {}
  ValueType GetValue() const override {
    return value_;
  }
};

template <class ValueType>
class ConditionFaceDerivative : public ConditionFace {
 public:
  virtual ValueType GetDerivative() const = 0;
};

template <class ValueType>
class ConditionFaceDerivativeFixed : public ConditionFaceDerivative<ValueType> {
  ValueType derivative_;
 public:
  explicit ConditionFaceDerivativeFixed(const ValueType& derivative)
      : derivative_(derivative)
  {}
  virtual ValueType GetDerivative() const override {
    return derivative_;
  }
};

class ConditionCell {
 public:
  virtual ~ConditionCell() {}
};

template <class ValueType>
class ConditionCellValue : public ConditionCell {
 public:
  virtual ValueType GetValue() const = 0;
};

template <class ValueType>
class ConditionCellValueFixed : public ConditionCellValue<ValueType> {
  ValueType value_;
 public:
  explicit ConditionCellValueFixed(const ValueType& value)
      : value_(value)
  {}
  ValueType GetValue() const override {
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
  const geom::FieldFace<Scal>& probe_;
  const Scal threshold_;
 public:
  InterpolationInnerFaceFirstUpwind(const Mesh& mesh,
                               const geom::FieldFace<Scal>& probe,
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
  const geom::FieldFace<Scal>& probe_;
  const geom::FieldCell<Scal>& fc_prev_;
  const geom::FieldCell<Vect>& fc_prev_grad_;
  const Scal threshold_;
 public:
  InterpolationInnerFaceSecondUpwindDeferred(
      const Mesh& mesh,
      const geom::FieldFace<Scal>& probe,
      const geom::FieldCell<Scal>& fc_prev,
      const geom::FieldCell<Vect>& fc_prev_grad,
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
  const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_cond_;
 public:
  InterpolationBoundaryFaceNearestCell(
      const Mesh& mesh,
      const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_cond)
      : Approximation<Mesh, IdxFace, Expr>(mesh)
      , mf_cond_(mf_cond)
  {}
  Expr GetExpression(IdxFace idxface) const override {
    const Mesh& mesh = this->mesh;
    Expr expr;
    if (auto cond_generic = mf_cond_.find(idxface)) {
      if (auto cond_value =
          dynamic_cast<ConditionFaceValue<Scal>*>(cond_generic->get())) {
        expr.SetConstant(cond_value->GetValue());
      } else if (auto cond_derivative =
          dynamic_cast<ConditionFaceDerivative<Scal>*>(cond_generic->get())) {
        size_t id = mesh.GetValidNeighbourCellId(idxface);
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
  const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_cond_;
 public:
  explicit DerivativeBoundaryFacePlain(
      const Mesh& mesh,
      const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_cond)
      : Approximation<Mesh, IdxFace, Expr>(mesh)
      , mf_cond_(mf_cond)
  {}
  Expr GetExpression(IdxFace idxface) const override {
    const Mesh& mesh = this->mesh;
    Expr expr;
    if (auto cond_generic = mf_cond_.find(idxface)) {
      if (auto cond_derivative =
          dynamic_cast<ConditionFaceDerivative<Scal>*>(cond_generic->get())) {
        expr.SetConstant(cond_derivative->GetDerivative());
      } else if (auto cond_value =
          dynamic_cast<ConditionFaceValue<Scal>*>(cond_generic->get())) {
        size_t id = mesh.GetValidNeighbourCellId(idxface);
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
geom::FieldFace<T> Interpolate(const geom::FieldNode<T>& fn_u,
                               const Mesh& mesh) {
  geom::FieldFace<T> res(mesh);

  for (auto idxface : mesh.Faces()) {
    T sum = static_cast<T>(0);
    for (size_t i = 0; i < mesh.GetNumNeighbourNodes(idxface); ++i) {
      sum += fn_u[mesh.GetNeighbourNode(idxface, i)];
    }
    res[idxface] = sum / mesh.GetNumNeighbourNodes(idxface);
  }

  return res;
}

template <class T, class Mesh>
T GetInterpolatedInner(const geom::FieldCell<T>& fc_u,
                       IdxFace idxface,
                       const Mesh& mesh) {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
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
geom::Vect<Scal, dim> GetGeometricAverage(const geom::Vect<Scal, dim>& a,
                                          const geom::Vect<Scal, dim>& b,
                                          Scal threshold = 1e-8) {
  auto arithmetic_aver = (a + b) * 0.5;
  if (arithmetic_aver.norm() < threshold) {
    return arithmetic_aver;
  }
  auto direction = arithmetic_aver / arithmetic_aver.norm();
  return direction * std::sqrt(a.norm() * b.norm());
}

template <class T, class Mesh>
geom::FieldFace<T> Interpolate(
    const geom::FieldCell<T>& fc_u,
    const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_cond_u,
    const Mesh& mesh, bool geometric = false) {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;

  geom::FieldFace<T> res(mesh, T(0)); // Valid value essential for extrapolation

  if (geometric) {
#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Faces().size()); ++i) {
      IdxFace idxface(i);
      if (mesh.IsInner(idxface)) {
        IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
        IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
        res[idxface] = GetGeometricAverage(fc_u[cm], fc_u[cp]);
      }
    }
  } else {
#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Faces().size()); ++i) {
      IdxFace idxface(i);
      if (mesh.IsInner(idxface)) {
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
  }

  for (auto it = mf_cond_u.cbegin(); it != mf_cond_u.cend(); ++it) {
    IdxFace idxface = it->GetIdx();
    if (mesh.IsInner(idxface)) {
      continue;
    }
    ConditionFace* cond = it->GetValue().get();
    if (auto cond_value = dynamic_cast<ConditionFaceValue<T>*>(cond)) {
      res[idxface] = cond_value->GetValue();
    } else if (auto cond_derivative =
        dynamic_cast<ConditionFaceDerivative<T>*>(cond)) {
      size_t id = mesh.GetValidNeighbourCellId(idxface);
      IdxCell cc = mesh.GetNeighbourCell(idxface, id);
      Scal factor = (id == 0 ? 1. : -1.);
      Scal alpha = mesh.GetVectToCell(idxface, id).norm() * factor;
      res[idxface] = fc_u[cc] + cond_derivative->GetDerivative() * alpha;
    } else if (dynamic_cast<ConditionFaceExtrapolation*>(cond)) {
      size_t id = mesh.GetValidNeighbourCellId(idxface);
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
geom::FieldCell<typename Mesh::Scal>
CalcLaplacian(const geom::FieldCell<typename Mesh::Scal>& fc_u,
              const Mesh& mesh) {
  using Scal = typename Mesh::Scal;
  geom::FieldCell<Scal> res(mesh);

  for (auto idxcell : mesh.Cells()) {
    Scal sum = 0.;
    for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
      auto idxface = mesh.GetNeighbourFace(idxcell, i);
      if (mesh.IsInner(idxface)) {
        auto cm = mesh.GetNeighbourCell(idxface, 0);
        auto cp = mesh.GetNeighbourCell(idxface, 1);
        auto dm = mesh.GetVectToCell(idxface, 0);
        auto dp = mesh.GetVectToCell(idxface, 1);
        Scal derivative = (fc_u[cp] - fc_u[cm]) / (dp - dm).norm();
        sum += derivative * mesh.GetArea(idxface) *
            mesh.GetOutwardFactor(idxcell, i);
      }
    }
    res[idxcell] = sum / mesh.GetVolume(idxcell);
  }

  return res;
}

template <class T, class Mesh>
geom::FieldFace<T> InterpolateFirstUpwind(
    const geom::FieldCell<T>& fc_u,
    const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_cond_u,
    const geom::FieldFace<typename Mesh::Scal>& probe,
    const Mesh& mesh, typename Mesh::Scal threshold = 1e-8) {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;

  geom::FieldFace<T> res(mesh);

  for (IdxFace idxface : mesh.Faces()) {
     if (mesh.IsInner(idxface)) {
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
  }

  for (auto it = mf_cond_u.cbegin(); it != mf_cond_u.cend(); ++it) {
    IdxFace idxface = it->GetIdx();
    if (mesh.IsInner(idxface)) {
      continue;
    }
    ConditionFace* cond = it->GetValue().get();
    if (auto cond_value = dynamic_cast<ConditionFaceValue<T>*>(cond)) {
      res[idxface] = cond_value->GetValue();
    } else if (auto cond_derivative =
        dynamic_cast<ConditionFaceDerivative<T>*>(cond)) {
      size_t id = mesh.GetValidNeighbourCellId(idxface);
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
geom::FieldFace<typename Mesh::Scal>
InterpolateSuperbee(
    const geom::FieldCell<typename Mesh::Scal>& fc_u,
    const geom::FieldCell<typename Mesh::Vect>& fc_u_grad,
    const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_cond_u,
    const geom::FieldFace<typename Mesh::Scal>& probe,
    const Mesh& mesh, typename Mesh::Scal threshold = 1e-8) {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;

  geom::FieldFace<Scal> res(mesh);

  for (IdxFace idxface : mesh.Faces()) {
     if (mesh.IsInner(idxface)) {
       IdxCell P = mesh.GetNeighbourCell(idxface, 0);
       IdxCell E = mesh.GetNeighbourCell(idxface, 1);
       Vect rp = mesh.GetVectToCell(idxface, 0); 
       Vect re = mesh.GetVectToCell(idxface, 1); 
       const auto& u = fc_u;
       const auto& g = fc_u_grad;
       if (probe[idxface] > threshold) {
         res[idxface] =
             u[P] + 0.5 * Superbee(u[E] - u[P],
                                   -4. * g[P].dot(rp) - (u[E] - u[P]));
       } else if (probe[idxface] < -threshold) {
         res[idxface] =
             u[E] - 0.5 * Superbee(u[E] - u[P],
                                   4. * g[E].dot(re) - (u[E] - u[P]));
       } else {
         // TODO: Make a CDS proper for non-uniform grids
         res[idxface] = 0.5 * (u[P] + u[E]);
       }
     }
  }

  // TODO: Move interpolation on boundaries to a function
  for (auto it = mf_cond_u.cbegin(); it != mf_cond_u.cend(); ++it) {
    IdxFace idxface = it->GetIdx();
    if (mesh.IsInner(idxface)) {
      continue;
    }
    ConditionFace* cond = it->GetValue().get();
    if (auto cond_value = dynamic_cast<ConditionFaceValue<Scal>*>(cond)) {
      res[idxface] = cond_value->GetValue();
    } else if (auto cond_derivative =
        dynamic_cast<ConditionFaceDerivative<Scal>*>(cond)) {
      size_t id = mesh.GetValidNeighbourCellId(idxface);
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
geom::FieldCell<T> Average(const geom::FieldFace<T>& ff_u, const Mesh& mesh) {
  using Scal = typename Mesh::Scal;
  geom::FieldCell<T> res(mesh);
  for (IdxCell idxcell : mesh.Cells()) {
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
geom::FieldCell<T>
GetSmoothField(const geom::FieldCell<T>& fc_u,
            const Mesh& mesh, size_t repeat = 1) {
  geom::MapFace<std::shared_ptr<ConditionFace>> mf_cond;

  for (auto idxface : mesh.Faces()) {
    if (!mesh.IsExcluded(idxface) && !mesh.IsInner(idxface)) {
      mf_cond[idxface] =
          std::make_shared<ConditionFaceDerivativeFixed<T>>(T(0));
    }
  }

  geom::FieldCell<T> res = fc_u;

  for (size_t i = 0; i < repeat; ++i) {
    res = Average(Interpolate(res, mf_cond, mesh), mesh);
  }

  return res;
}

template <class Mesh>
geom::FieldCell<typename Mesh::Vect> Gradient(
    const geom::FieldFace<typename Mesh::Scal>& ff_u,
    const Mesh& mesh) {
  using Vect = typename Mesh::Vect;
  geom::FieldCell<Vect> res(mesh, Vect::kZero);
#pragma omp parallel for
  for (IntIdx rawcell = 0; rawcell < static_cast<IntIdx>(mesh.Cells().size()); ++rawcell) {
    IdxCell idxcell(rawcell);
    if (!mesh.IsExcluded(idxcell)) {
      Vect sum = Vect::kZero;
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
        sum += mesh.GetOutwardSurface(idxcell, i) * ff_u[idxface];
      }
      res[idxcell] = sum / mesh.GetVolume(idxcell);
    }
  }
  return res;
}

class UnsteadySolver {
  double time_;
  double time_step_;
 public:
  UnsteadySolver(double time, double time_step)
      : time_(time)
      , time_step_(time_step)
  {}
  virtual ~UnsteadySolver() {}
  virtual double GetTime() const {
    return time_;
  }
  virtual void SetTime(double time) {
    time_  = time;
  }
  virtual double GetTimeStep() const {
    return time_step_;
  }
  virtual void SetTimeStep(double time_step) {
    time_step_ = time_step;
  }
  virtual void IncTime() {
    time_ += time_step_;
  }
  virtual void StartStep() {}
  virtual void CalcStep() = 0;
  virtual void FinishStep() {
    IncTime();
  }
};

class UnsteadyIterativeSolver : public UnsteadySolver {
  size_t iteration_count_;
  double convergence_tolerance_;
  size_t num_iterations_limit_;

 protected:
  void IncIterationCount() {
    ++iteration_count_;
  }
  void ClearIterationCount() {
    iteration_count_ = 0;
  }

 public:
  UnsteadyIterativeSolver(double time, double time_step,
                          double convergence_tolerance,
                          size_t num_iterations_limit)
      : UnsteadySolver(time, time_step)
      , iteration_count_(0)
      , convergence_tolerance_(convergence_tolerance)
      , num_iterations_limit_(num_iterations_limit)
  {}
  virtual void MakeIteration() = 0;
  virtual bool IsConverged() const {
    return this->GetIterationCount() >= num_iterations_limit_ ||
        GetConvergenceIndicator() < convergence_tolerance_;
  }
  virtual double GetConvergenceIndicator() const {
    return 0.;
  }
  virtual size_t GetIterationCount() const {
    return iteration_count_;
  }
  virtual void StartStep() override {
    ClearIterationCount();
  }
  virtual void CalcStep() override {
    while (!IsConverged()) {
      MakeIteration();
    }
  }
};

// After FinishStep() call:
//   iter_curr, iter_prev -- undefined
//   time_curr -- latest state (n)
//   time_prev -- previous state (n-1)
// Between StartStep() and FinishStep():
//   iter_curr -- latest iteration having been calculated (n+1, s)
//   iter_prev -- previous iteration (n+1, s-1)
//   time_curr -- preceding the state being calculated (n)
//   time_prev -- previous state (n-1)
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
  using Idx = typename Field::IdxType;
  for (Idx idx : geom::Range<Idx>(mesh)) {
    res = std::max(res, std::abs(first[idx] - second[idx]));
  }
  return res;
}

template <class Idx, class Mesh, class Scal = typename Mesh::Scal>
Scal CalcDiff(const geom::FieldGeneric<typename Mesh::Vect, Idx>& first,
                const geom::FieldGeneric<typename Mesh::Vect, Idx>& second,
                const Mesh& mesh) {
  Scal res = 0.;
  for (Idx idx : geom::Range<Idx>(mesh)) {
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
