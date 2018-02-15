/*
 * advection.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "mesh.hpp"
#include "linear.hpp"
#include <exception>
#include "solver.hpp"
#include "particle_system.hpp"

#include <fstream>

namespace solver {

template <class Mesh>
geom::FieldCell<typename Mesh::Vect> GetDeformingVelocity(const Mesh& mesh) {
  using Vect = typename Mesh::Vect;
  geom::FieldCell<Vect> res(mesh, Vect::kZero);
  for (auto idxcell : mesh.Cells()) {
    auto x = mesh.GetCenter(idxcell);
    res[idxcell][0] = -std::cos(x[1]) * std::sin(x[0]);
    res[idxcell][1] = std::cos(x[0]) * std::sin(x[1]);
  }
  return res;
}

template <class Mesh, class VelocityField>
class AdvectionSolver : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;

 protected:
  const VelocityField* p_f_velocity_;
  const geom::FieldCell<Scal>* p_fc_source_;

 public:
  AdvectionSolver(double time, double time_step,
                  const VelocityField* p_f_velocity,
                  const geom::FieldCell<Scal>* p_fc_source)
      : UnsteadyIterativeSolver(time, time_step, 1e-5, 1)
      , p_f_velocity_(p_f_velocity)
      , p_fc_source_(p_fc_source)
  {}
  virtual void AssignVelocity(const VelocityField* p_f_velocity) {
    p_f_velocity_ = p_f_velocity;
  }
  virtual const geom::FieldCell<Scal>& GetField() = 0;
};
/*
template <class Mesh, class VelocityField>
class AdvectionFactory {
 public:
  virtual std::shared_ptr<AdvectionSolver<Mesh, VelocityField>> Create() = 0;
  virtual ~AdvectionFactory() {}
};
*/
template <class Mesh, class VelocityField>
class AdvectionSolverMulti : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;

 protected:
  const VelocityField* p_f_velocity_;
  std::vector<const geom::FieldCell<Scal>*> v_p_fc_source_;

 public:
  AdvectionSolverMulti(double time, double time_step,
                       const VelocityField* p_f_velocity,
                       const std::vector<const geom::FieldCell<Scal>*>&
                       v_p_fc_source)
      : UnsteadyIterativeSolver(time, time_step, 1e-5, 1)
      , p_f_velocity_(p_f_velocity)
      , v_p_fc_source_(v_p_fc_source)
  {}
  virtual void AssignVelocity(const VelocityField* p_f_velocity) {
    p_f_velocity_ = p_f_velocity;
  }
  virtual const geom::FieldCell<Scal>& GetField(size_t field_id) = 0;
};
/*
template <class Mesh, class VelocityField>
class AdvectionSolverMultiFromSingle :
    AdvectionSolverMulti<Mesh, VelocityField> {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  std::vector<std::shared_ptr<AdvectionSolver<Mesh, VelocityField>>> solvers_;
 public:
  AdvectionSolverMultiFromSingle(
      size_t num_fields,
      const VelocityField* p_f_velocity,
      AdvectionFactory<Mesh, VelocityField>& factory,
      double time, double time_step)
      : AdvectionSolverMulti<Mesh, VelocityField>(
          num_fields, p_f_velocity, time, time_step)
  {}
  void AssignVelocity(const VelocityField* p_f_velocity) {
    for (auto& solver : solvers_) {
      solver->AssignVelocity(p_f_velocity);
    }
  }
  const geom::FieldCell<Scal>& GetField(size_t field_id) {
    return solvers_[field_id]->GetField();
  }
};
*/
template <class Mesh, class VelocityField>
class AdvectionSolverImplicit :
    public AdvectionSolver<Mesh, VelocityField> {
  const Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  LayersData<geom::FieldCell<Scal>> fc_u_;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  geom::MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;
  std::shared_ptr<LinearSolver<Scal, IdxCell, Expr>> linear_;
  // Common buffers:
  geom::FieldFace<Vect> ff_velocity_;
  geom::FieldFace<Scal> ff_volume_flux_;
  geom::FieldFace<Expr> ff_flux_;
  geom::FieldCell<Expr> fc_flux_sum_, fc_system_;

  bool time_second_order_;

 public:
  AdvectionSolverImplicit(
      const Mesh& mesh,
      const geom::FieldCell<Scal>& fc_u_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond,
      const VelocityField* p_fn_velocity,
      const geom::FieldCell<Scal>* p_fc_source,
      double time, double time_step,
      const LinearSolverFactory& linear_factory,
      bool time_second_order = false)
      : AdvectionSolver<Mesh, VelocityField>(
          time, time_step, p_fn_velocity, p_fc_source)
      , mesh(mesh)
      , mf_u_cond_(mf_u_cond)
      , ff_volume_flux_(mesh)
      , ff_flux_(mesh)
      , fc_flux_sum_(mesh)
      , fc_system_(mesh)
      , time_second_order_(time_second_order)
  {
    fc_u_.time_curr = fc_u_initial;
    fc_u_.time_prev = fc_u_initial;

    linear_ = linear_factory.Create<Scal, IdxCell, Expr>();

    for (auto it = mf_u_cond_.cbegin() ;
        it != mf_u_cond_.cend();
        ++it) {
        mf_u_cond_[it->GetIdx()] = it->GetValue().get();
    }
  }
  void StartStep() override {
    this->ClearIterationCount();
    fc_u_.iter_curr = fc_u_.time_curr;
  }
  geom::FieldFace<Scal> ConvertVolumeFlux(
      const geom::FieldFace<Scal>* p_f_velocity) {
    return *p_f_velocity;
  }
  geom::FieldFace<Scal> ConvertVolumeFlux(
      const geom::FieldNode<Vect>* p_f_velocity)  {
    auto ff_velocity = Interpolate(*p_f_velocity, mesh);

    geom::FieldFace<Scal> ff_volume_flux(mesh);
    for (IdxFace idxface : mesh.Faces()) {
      ff_volume_flux[idxface] =
          ff_velocity[idxface].dot(mesh.GetSurface(idxface));
    }
    return ff_volume_flux;
  }
  void MakeIteration() override {
    auto& prev = fc_u_.iter_prev;
    auto& curr = fc_u_.iter_curr;
    prev = curr;

    ff_volume_flux_ = ConvertVolumeFlux(this->p_f_velocity_);

    InterpolationInnerFaceFirstUpwind<Mesh, Expr>
    value_inner(mesh, ff_volume_flux_);

    InterpolationBoundaryFaceNearestCell<Mesh, Expr>
    value_boundary(mesh, mf_u_cond_);

    for (auto idxface : mesh.Faces()) {
      if (mesh.IsInner(idxface)) {
        ff_flux_[idxface] =
            value_inner.GetExpression(idxface) * ff_volume_flux_[idxface];
      } else {
        ff_flux_[idxface] =
            value_boundary.GetExpression(idxface) * ff_volume_flux_[idxface];
      }
    }

    for (auto idxcell : mesh.Cells()) {
      auto& expr = fc_flux_sum_[idxcell];
      expr.Clear();
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
        expr += ff_flux_[idxface] * mesh.GetOutwardFactor(idxcell, i);
      }
    }

    bool implicit = true;

    auto dt = this->GetTimeStep();
    auto coeffs = GetDerivativeApproxCoeffs(
        0., {-2. * dt, -dt, 0.}, time_second_order_ ? 0 : 1);

    for (auto idxcell : mesh.Cells()) {
      const auto& flux = fc_flux_sum_[idxcell];
      auto vol = mesh.GetVolume(idxcell);
      Scal diag = vol * coeffs[2];
      Expr unsteady;
      unsteady.InsertTerm(diag, idxcell);
      unsteady.SetConstant(
          vol * coeffs[0] * fc_u_.time_prev[idxcell] +
          vol * coeffs[1] * fc_u_.time_curr[idxcell] +
          -mesh.GetVolume(idxcell) * (*this->p_fc_source_)[idxcell]);
      auto& expr = fc_system_[idxcell];
      expr = unsteady + (implicit ? flux : Expr(flux.Evaluate(prev)));
    }

    for (auto idxcell : mesh.Cells()) {
      auto& expr = fc_system_[idxcell];
      expr.SetConstant(expr.Evaluate(prev));
    }

    curr = linear_->Solve(fc_system_);
    for (auto idxcell : mesh.Cells()) {
      curr[idxcell] += prev[idxcell];
    }

    this->IncIterationCount();
  }
  void FinishStep() override {
    fc_u_.time_prev = fc_u_.time_curr;
    fc_u_.time_curr = fc_u_.iter_curr;
    this->IncTime();
  }
  double GetConvergenceIndicator() const {
    if (this->GetIterationCount() == 0) {
      return 1.;
    }
    return CalcDiff<geom::FieldCell<Scal>, Mesh>(
        fc_u_.iter_curr, fc_u_.iter_prev, mesh);
  }
  const geom::FieldCell<Scal>& GetField() override {
    return fc_u_.time_curr;
  }
};

template <class Mesh, class VelocityField>
class AdvectionSolverExplicit :
    public AdvectionSolver<Mesh, VelocityField> {
  const Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  LayersData<geom::FieldCell<Scal>> fc_u_;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  geom::MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;
  // Common buffers:
  geom::FieldFace<Vect> ff_velocity_;
  geom::FieldFace<Scal> ff_volume_flux_;
  geom::FieldFace<Scal> ff_flux_;
  geom::FieldFace<Scal> ff_u_;

 public:
  AdvectionSolverExplicit(
      const Mesh& mesh,
      const geom::FieldCell<Scal>& fc_u_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond_,
      const VelocityField* p_fn_velocity,
      const geom::FieldCell<Scal>* p_fc_source,
      double time, double time_step)
      : AdvectionSolver<Mesh, VelocityField>(
          time, time_step, p_fn_velocity, p_fc_source)
      , mesh(mesh)
      , mf_u_cond_(mf_u_cond_)
      , ff_volume_flux_(mesh)
      , ff_flux_(mesh)
  {
    fc_u_.time_curr = fc_u_initial;
  }
  void StartStep() override {
    this->ClearIterationCount();
    fc_u_.time_prev = fc_u_.time_curr;
    fc_u_.iter_curr = fc_u_.time_prev;
  }
  // Correct the inconsistency with arguments: velocity vs volume flux
  geom::FieldFace<Scal> ConvertVolumeFlux(
      const geom::FieldFace<Scal>* p_f_velocity) {
    return *p_f_velocity;
  }
  geom::FieldFace<Scal> ConvertVolumeFlux(
      const geom::FieldNode<Vect>* p_f_velocity)  {
    auto ff_velocity = Interpolate(*p_f_velocity, mesh);

    geom::FieldFace<Scal> ff_volume_flux(mesh);
    for (IdxFace idxface : mesh.Faces()) {
      ff_volume_flux[idxface] =
          ff_velocity[idxface].dot(mesh.GetSurface(idxface));
    }
    return ff_volume_flux;
  }
  void MakeIteration() override {
    auto& prev = fc_u_.iter_prev;
    auto& curr = fc_u_.iter_curr;
    prev = curr;

    ff_volume_flux_ = ConvertVolumeFlux(this->p_f_velocity_);

    ff_u_ = InterpolateSuperbee(
        prev,
        Gradient(Interpolate(prev, mf_u_cond_, mesh), mesh),
        mf_u_cond_, ff_volume_flux_, mesh);

    for (auto idxface : mesh.Faces()) {
      ff_flux_[idxface] = ff_u_[idxface] * ff_volume_flux_[idxface];
    }

    for (auto idxcell : mesh.Cells()) {
      Scal flux_sum = 0.;
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
        flux_sum += ff_flux_[idxface] * mesh.GetOutwardFactor(idxcell, i);
      }

      curr[idxcell] = fc_u_.time_prev[idxcell] -
          this->GetTimeStep() / mesh.GetVolume(idxcell) * flux_sum +
          this->GetTimeStep() * (*this->p_fc_source_)[idxcell];
    }

    this->IncIterationCount();
  }
  void FinishStep() override {
    fc_u_.time_curr = fc_u_.iter_curr;
    this->IncTime();
  }
  double GetConvergenceIndicator() const override {
    if (this->GetIterationCount() == 0) {
      return 1.;
    }
    return CalcDiff<geom::FieldCell<Scal>, Mesh>(
        fc_u_.iter_curr, fc_u_.iter_prev, mesh);
  }
  const geom::FieldCell<Scal>& GetField() override {
    return fc_u_.time_curr;
  }
};


template <class Mesh, class VelocityField>
class AdvectionSolverMultiExplicit :
    public AdvectionSolverMulti<Mesh, VelocityField> {
  const Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  size_t num_fields_;
  std::vector<LayersData<geom::FieldCell<Scal>>> v_fc_u_;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  std::vector<geom::MapFace<std::shared_ptr<ConditionFace>>> v_mf_u_cond_;
  std::vector<const VelocityField*> v_p_ff_volume_flux_slip_;
  bool spatial_split_;
  Scal sharp_;
  std::vector<Scal> sharp_max_;
  // Common buffers:
  geom::FieldFace<Scal> ff_u_;

 public:
  AdvectionSolverMultiExplicit(
      const Mesh& mesh,
      const std::vector<geom::FieldCell<Scal>>& v_fc_u_initial,
      const std::vector<geom::MapFace<std::shared_ptr<ConditionFace>>>&
      v_mf_u_cond,
      const VelocityField* p_fn_velocity,
      const std::vector<const VelocityField*>& v_p_ff_volume_flux_slip,
      const std::vector<const geom::FieldCell<Scal>*>& v_p_fc_source,
      double time, double time_step,
      bool spatial_split, Scal sharp, std::vector<Scal> sharp_max)
      : AdvectionSolverMulti<Mesh, VelocityField>(
          time, time_step, p_fn_velocity, v_p_fc_source)
      , mesh(mesh)
      , num_fields_(v_fc_u_initial.size())
      , v_mf_u_cond_(v_mf_u_cond)
      , v_p_ff_volume_flux_slip_(v_p_ff_volume_flux_slip)
      , spatial_split_(spatial_split)
      , sharp_(sharp)
      , sharp_max_(sharp_max)
  {
    v_fc_u_.resize(num_fields_);
    v_mf_u_cond_.resize(num_fields_);
    for (size_t field = 0; field < num_fields_; ++field) {
      v_fc_u_[field].time_curr = v_fc_u_initial[field];
    }

  }
  void StartStep() override {
    this->ClearIterationCount();
    for (size_t field = 0; field < num_fields_; ++field) {
      v_fc_u_[field].time_prev = v_fc_u_[field].time_curr;
      v_fc_u_[field].iter_curr = v_fc_u_[field].time_prev;
    }
  }
  // Correct the inconsistency with arguments: velocity vs volume flux
  geom::FieldFace<Scal> ConvertVolumeFlux(
      const geom::FieldFace<Scal>* p_f_velocity) {
    return *p_f_velocity;
  }
  geom::FieldFace<Scal> ConvertVolumeFlux(
      const geom::FieldNode<Vect>* p_f_velocity)  {
    auto ff_velocity = Interpolate(*p_f_velocity, mesh);

    geom::FieldFace<Scal> ff_volume_flux(mesh);
    for (IdxFace idxface : mesh.Faces()) {
      ff_volume_flux[idxface] =
          ff_velocity[idxface].dot(mesh.GetSurface(idxface));
    }
    return ff_volume_flux;
  }
  void MakeIteration() override {
    const auto& ff_volume_flux_mixture_ =
        ConvertVolumeFlux(this->p_f_velocity_);

    for (size_t field = 0; field < num_fields_; ++field) {
      auto& prev = v_fc_u_[field].iter_prev;
      auto& curr = v_fc_u_[field].iter_curr;
      prev = curr;

      const auto& ff_volume_flux_slip =
          ConvertVolumeFlux(v_p_ff_volume_flux_slip_[field]);
      auto ff_volume_flux = ff_volume_flux_mixture_;
      for (auto idxface : mesh.Faces()) {
        ff_volume_flux[idxface] += ff_volume_flux_slip[idxface];
      }

      // Apply operator in each direction
      size_t num_stages = (spatial_split_ ? dim : 1);
      for (size_t stage = 0; stage < num_stages; ++stage) {
        ff_u_ = InterpolateSuperbee(
            curr,
            Gradient(Interpolate(curr, v_mf_u_cond_[field], mesh), mesh),
            v_mf_u_cond_[field], ff_volume_flux, mesh);

//        ff_u_ = solver::InterpolateFirstUpwind(
//            curr, v_mf_u_cond_[field], ff_volume_flux, mesh);
        for (auto idxcell : mesh.Cells()) {
          Scal flux_sum = 0.;
          for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
            IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
            if ((i / 2) % num_stages == stage) {
              flux_sum += ff_u_[idxface] * ff_volume_flux[idxface] *
                  mesh.GetOutwardFactor(idxcell, i);
            }
          }

          curr[idxcell] +=
              -this->GetTimeStep() / mesh.GetVolume(idxcell) * flux_sum;
        }
      }

      // Interface sharpening
      // zero-derivative bc for Vect
      geom::MapFace<std::shared_ptr<ConditionFace>> mfvz;
      for (auto idxface : mesh.Faces()) {
        if (!mesh.IsExcluded(idxface) && !mesh.IsInner(idxface)) {
          mfvz[idxface] =
              std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect(0));
        }
      }
      auto af = Interpolate(curr, v_mf_u_cond_[field], mesh);
      auto gc = Gradient(af, mesh);
      auto gf = Interpolate(gc, mfvz, mesh);
      geom::FieldFace<Scal> ff(mesh, 0);
      if (std::abs(sharp_) > 1e-10) {
        for (auto idxface : mesh.Faces()) {
          auto n = gf[idxface];
          n /= (n.norm() + 1e-6);
          auto nf = n.dot(mesh.GetNormal(idxface));
          auto uf = ff_volume_flux[idxface];
          auto am = sharp_max_[field];
          auto epsh = sharp_ * mesh.GetArea(idxface);
          ff[idxface] = std::abs(uf * nf) * nf * 
            (epsh * gf[idxface].norm() - af[idxface] * 
                  (1. - af[idxface] / am));
        }
      }

      // Apply sources
      for (auto idxcell : mesh.Cells()) {
        curr[idxcell] +=
            this->GetTimeStep() * (*this->v_p_fc_source_[field])[idxcell];

        // sharpening
        Scal sh = 0.;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          sh += mesh.GetOutwardFactor(idxcell, i) * ff[idxface];
        }
        curr[idxcell] += this->GetTimeStep() * sh / mesh.GetVolume(idxcell);
      }

      // Normalize
      //for (auto idxcell : mesh.Cells()) {
      //  auto& c = curr[idxcell];
      //  c = std::max(0., std::min(sharp_max_[field], c));
      //}
    }
    this->IncIterationCount();
  }
  void FinishStep() override {
    for (size_t field = 0; field < num_fields_; ++field) {
      v_fc_u_[field].time_curr = v_fc_u_[field].iter_curr;
    }
    this->IncTime();
  }
  double GetConvergenceIndicator() const override {
    if (this->GetIterationCount() == 0) {
      return 1.;
    }
    return 0.;
  }
  const geom::FieldCell<Scal>& GetField(size_t field_id) override {
    return v_fc_u_[field_id].time_curr;
  }
};
/*
template <class Mesh>
class AdvectionSolverParticles :
    public AdvectionSolver<Mesh, geom::FieldNode<typename Mesh::Vect>> {
  const Mesh& mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using VelocityField = geom::FieldNode<typename Mesh::Vect>;
  geom::MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;
  geom::MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;
  ParticleSystem<Mesh> particles_;

 public:
  AdvectionSolverParticles(
      const Mesh& mesh,
      const geom::FieldCell<Scal>& fc_u_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond,
      geom::FieldNode<Vect>* p_fn_velocity,
      double time, double time_step)
      : AdvectionSolver<Mesh, VelocityField>(time, time_step, p_fn_velocity)
      , mesh(mesh)
      , mf_u_cond_(mf_u_cond)
      , particles_(mesh, {fc_u_initial}, p_fn_velocity, time, time_step)
  {
    for (auto it = mf_u_cond_.cbegin() ;
        it != mf_u_cond_.cend();
        ++it) {
        mf_u_cond_[it->GetIdx()] = it->GetValue().get();
    }
  }
  void StartStep() override {
    this->ClearIterationCount();
  }
  void MakeIteration() override {
    particles_.CalcStep();

    this->IncIterationCount();
  }
  void FinishStep() override {
    this->IncTime();
  }
  double GetConvergenceIndicator() const {
    return 1.;
  }
  const geom::FieldCell<Scal>& GetField() override {
    return particles_.GetField(0);
  }
};
*/

template <class Mesh, class VelocityField>
class AdvectionSolverMultiParticles :
    public AdvectionSolverMulti<Mesh, VelocityField> {
  static constexpr size_t dim = Mesh::dim;
  const Mesh& mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  size_t num_fields_;
  using MIdx = geom::MIdxGeneral<dim>;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  std::vector<geom::MapFace<std::shared_ptr<ConditionFace>>>
  v_mf_u_cond_;
  geom::FieldNode<Vect> fn_velocity_;
  std::shared_ptr<ParticleSystem<Mesh>> particles_;
  std::vector<geom::FieldCell<Scal>> v_fc_source_extra_;
  // Common buffers:
  geom::FieldFace<Vect> ff_velocity_;
  geom::FieldFace<Scal> ff_volume_flux_;
  geom::FieldFace<Scal> ff_flux_;

 public:
  AdvectionSolverMultiParticles(
      const Mesh& mesh,
      const std::vector<geom::FieldCell<Scal>>& v_fc_u_initial,
      const std::vector<geom::MapFace<std::shared_ptr<ConditionFace>>>&
      v_mf_u_cond,
      const VelocityField* p_f_velocity,
      const std::vector<const geom::FieldCell<Scal>*>& v_p_fc_source,
      double time, double time_step,
      Scal spawning_gap_relative, Scal particle_radius_factor,
      size_t min_num_particle_in_cell, size_t max_num_particle_in_cell,
      Scal back_relaxation_factor)
      : AdvectionSolverMulti<Mesh, VelocityField>(
          time, time_step, p_f_velocity, v_p_fc_source)
      , mesh(mesh)
      , num_fields_(v_fc_u_initial.size())
      , v_mf_u_cond_(v_mf_u_cond)
      , v_fc_source_extra_(num_fields_)
  {

    std::vector<const geom::FieldCell<Scal>*> v_p_fc_source_extra(num_fields_);
    for (size_t i = 0; i < num_fields_; ++i) {
      v_fc_source_extra_[i].Reinit(mesh);
      v_p_fc_source_extra[i] = &v_fc_source_extra_[i];
    }

    particles_ = std::make_shared<ParticleSystem<Mesh>>(
        mesh, v_fc_u_initial, v_p_fc_source_extra, &fn_velocity_,
        time, time_step,
        spawning_gap_relative, particle_radius_factor,
        min_num_particle_in_cell, max_num_particle_in_cell,
        back_relaxation_factor);
  }
  void StartStep() override {
    this->ClearIterationCount();
  }
  geom::FieldNode<Vect> GetVelocity(
      const geom::FieldFace<Scal>* p_f_velocity) {
    geom::FieldCell<Vect> fc_velocity(mesh);
    for (auto idxcell : mesh.Cells()) {
      Vect sum = Vect::kZero;
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        auto idxface = mesh.GetNeighbourFace(idxcell, i);
        sum += mesh.GetCenter(idxface) * (*p_f_velocity)[idxface] *
            mesh.GetOutwardFactor(idxcell, i);
      }
      fc_velocity[idxcell] = sum / mesh.GetVolume(idxcell);
    }

    geom::FieldNode<Vect> res(mesh, Vect::kZero);
    geom::FieldNode<Scal> denominator(mesh, 0.);

    for (auto idxcell : mesh.Cells()) {
      for (size_t i = 0; i < mesh.GetNumNeighbourNodes(idxcell); ++i) {
        auto idxnode = mesh.GetNeighbourNode(idxcell, i);
        res[idxnode] += fc_velocity[idxcell];
        denominator[idxnode] += 1.;
      }
    }

    for (auto idxnode : mesh.Nodes()) {
      res[idxnode] /= denominator[idxnode];
    }

    return res;
  }
  geom::FieldNode<Vect> GetVelocityCartesianMesh(
      const geom::FieldFace<Scal>* p_f_velocity) {
    geom::FieldNode<Vect> res(mesh, Vect::kZero);
    geom::FieldNode<Vect> denominator(mesh, Vect::kZero);

    for (auto idxface : mesh.Faces()) {
      for (size_t i = 0; i < mesh.GetNumNeighbourNodes(idxface); ++i) {
        auto idxnode = mesh.GetNeighbourNode(idxface, i);
        auto dir = mesh.GetDirection(idxface);
        res[idxnode][dir] += (*p_f_velocity)[idxface] / mesh.GetArea(idxface);
        denominator[idxnode][dir] += 1.;
      }
    }

    for (auto idxnode : mesh.Nodes()) {
      res[idxnode][0] /= denominator[idxnode][0];
      res[idxnode][1] /= denominator[idxnode][1];
    }

    return res;
  }
  geom::FieldNode<Vect> GetVelocity(
      const geom::FieldNode<Vect>* p_f_velocity)  {
    return *p_f_velocity;
  }
  void MakeIteration() override {
    fn_velocity_ = GetVelocityCartesianMesh(this->p_f_velocity_);

    for (auto idxcell : mesh.Cells()) {
      Scal div = 0.;
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        auto idxface = mesh.GetNeighbourFace(idxcell, i);
        div += (*this->p_f_velocity_)[idxface] *
            mesh.GetOutwardFactor(idxcell, i);
      }
      div /= mesh.GetVolume(idxcell);

      for (size_t n = 0; n < num_fields_; ++n) {
        v_fc_source_extra_[n][idxcell] =
            (*this->v_p_fc_source_[n])[idxcell] - div * GetField(n)[idxcell];
      }
    }

    particles_->CalcStep();

    this->IncIterationCount();
  }
  void FinishStep() override {
    this->IncTime();
  }
  double GetConvergenceIndicator() const override {
    return 1.;
  }
  const geom::FieldCell<Scal>& GetField(size_t field_id) override {
    return particles_->GetField(field_id);
  }
};


} // namespace solver

