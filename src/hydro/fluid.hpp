#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"
#include "metrics.hpp"
#include "conv_diff.hpp"

namespace solver {

// TODO: Pass parameters as constructor arguments

template <class Mesh>
class ConvectionDiffusion : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  geom::FieldCell<Scal>* p_fc_density_;
  geom::FieldFace<Scal>* p_ff_kinematic_viscosity_; // adhoc: dynamic viscosity
  geom::FieldCell<Vect>* p_fc_force_;
  geom::FieldFace<Scal>* p_ff_vol_flux_;

 public:
  ConvectionDiffusion(double time, double time_step,
                      geom::FieldCell<Scal>* p_fc_density,
                      geom::FieldFace<Scal>* p_ff_kinematic_viscosity,
                      geom::FieldCell<Vect>* p_fc_force,
                      geom::FieldFace<Scal>* p_ff_vol_flux,
                      double convergence_tolerance,
                      size_t num_iterations_limit)
      : UnsteadyIterativeSolver(time, time_step, convergence_tolerance,
                                num_iterations_limit)
      , p_fc_density_(p_fc_density)
      , p_ff_kinematic_viscosity_(p_ff_kinematic_viscosity)
      , p_fc_force_(p_fc_force)
      , p_ff_vol_flux_(p_ff_vol_flux)
  {

  }
  virtual const geom::FieldCell<Vect>& GetVelocity() = 0;
  virtual const geom::FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual void CorrectVelocity(Layers layer,
                               const geom::FieldCell<Vect>& fc_corr) = 0;
  virtual const geom::FieldCell<Expr>& GetVelocityEquations(size_t comp) = 0;
};

template <class Mesh>
class ConvectionDiffusionImplicit : public ConvectionDiffusion<Mesh> {
  Mesh& mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using Solver = ConvectionDiffusionScalarImplicit<Mesh>;

  static constexpr size_t dim = Mesh::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  template <class T>
  using VectGeneric = std::array<T, dim>;
  LayersData<geom::FieldCell<Vect>> fc_velocity_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_velocity_cond_;
  geom::MapCell<std::shared_ptr<ConditionCell>> mc_velocity_cond_;
  geom::FieldCell<Vect>* p_fc_force_;

  VectGeneric<geom::MapFace<std::shared_ptr<ConditionFace>>>
  v_mf_velocity_cond_;
  // TODO: Extract scalar CellCondition
  VectGeneric<std::shared_ptr<Solver>>
  v_solver_;
  VectGeneric<geom::FieldCell<Scal>>
  v_fc_force_;

 public:
  void CopyToVector(Layers layer) {
    fc_velocity_.Get(layer).Reinit(mesh);
    for (size_t n = 0; n < dim; ++n) {
      SetComponent(fc_velocity_.Get(layer), n, v_solver_[n]->GetField(layer));
    }
  }
  ConvectionDiffusionImplicit(
      Mesh& mesh,
      const geom::FieldCell<Vect>& fc_velocity_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>&
      mf_velocity_cond,
      const geom::MapCell<std::shared_ptr<ConditionCell>>&
      mc_velocity_cond,
      Scal relaxation_factor,
      geom::FieldCell<Scal>* p_fc_density,
      geom::FieldFace<Scal>* p_ff_kinematic_viscosity,
      geom::FieldCell<Vect>* p_fc_force,
    geom::FieldFace<Scal>* p_ff_vol_flux,
    double time, double time_step,
    const LinearSolverFactory& linear_factory,
    double convergence_tolerance,
    size_t num_iterations_limit,
    bool time_second_order = true,
    Scal guess_extrapolation = 0.)
    : ConvectionDiffusion<Mesh>(
        time, time_step, p_fc_density, p_ff_kinematic_viscosity, p_fc_force,
        p_ff_vol_flux,
          convergence_tolerance, num_iterations_limit)
      , mesh(mesh)
      , mf_velocity_cond_(mf_velocity_cond)
      , mc_velocity_cond_(mc_velocity_cond)
      , p_fc_force_(p_fc_force)
  {
    for (size_t n = 0; n < dim; ++n) {
      // Boundary conditions for each velocity component
      // (copied from given vector conditions)
      for (auto it = mf_velocity_cond_.cbegin();
          it != mf_velocity_cond_.cend(); ++it) {
        IdxFace idxface = it->GetIdx();
        if (auto cond = dynamic_cast<ConditionFaceValue<Vect>*>(
            mf_velocity_cond_[idxface].get())) {
          v_mf_velocity_cond_[n][idxface] =
              std::make_shared<ConditionFaceValueExtractComponent<Vect>>(
                  cond, n);
        } else {
          throw std::runtime_error("Unknown boudnary condition type");
        }
      }

      // Initialize solver
      v_solver_[n] = std::make_shared<Solver>(
          mesh, GetComponent(fc_velocity_initial, n),
          v_mf_velocity_cond_[n],
          geom::MapCell<std::shared_ptr<ConditionCell>>() /*empty*/,
          relaxation_factor, p_fc_density, p_ff_kinematic_viscosity,
          &(v_fc_force_[n]), p_ff_vol_flux, time, time_step,
          linear_factory, convergence_tolerance,
          num_iterations_limit, time_second_order, guess_extrapolation);
    }
    CopyToVector(Layers::time_curr);
    CopyToVector(Layers::time_prev);
  }
  void StartStep() override {
    auto sem = mesh.GetSem("convdiffmulti-start");
    for (size_t n = 0; n < dim; ++n) {
      if (sem("dir-init")) {
        v_solver_[n]->SetTimeStep(this->GetTimeStep());
      }
      if (sem.Nested("dir-start")) {
        v_solver_[n]->StartStep();
      }
    }
    if (sem("tovect")) {
      CopyToVector(Layers::iter_curr);
      this->ClearIterationCount();
    }
  }
  void MakeIteration() override {
    auto sem = mesh.GetSem("convdiffmulti-iter");
    for (size_t n = 0; n < dim; ++n) {
      if (sem("dir-get")) {
        v_fc_force_[n] = GetComponent(*p_fc_force_, n);
      }
    }

    for (size_t n = 0; n < dim; ++n) {
      if (sem.Nested("dir-iter")) {
        v_solver_[n]->MakeIteration();
      }
    }

    if (sem("tovect")) {
      CopyToVector(Layers::iter_prev);
      CopyToVector(Layers::iter_curr);
      this->IncIterationCount();
    }
  }
  void FinishStep() override {
    auto sem = mesh.GetSem("convdiffmulti-finish");

    for (size_t n = 0; n < dim; ++n) {
      if (sem.Nested("dir-finish")) {
        v_solver_[n]->FinishStep();
      }
    }
    if (sem("tovect")) {
      CopyToVector(Layers::time_prev);
      CopyToVector(Layers::time_curr);
      this->IncTime();
    }
  }
  double GetConvergenceIndicator() const override {
    if (this->GetIterationCount() == 0) {
      return 1.;
    }
    return CalcDiff(fc_velocity_.iter_curr, fc_velocity_.iter_prev, mesh);
  }
  const geom::FieldCell<Vect>& GetVelocity() override {
    return fc_velocity_.time_curr;
  }
  const geom::FieldCell<Vect>& GetVelocity(Layers layer) override {
    return fc_velocity_.Get(layer);
  }
  void CorrectVelocity(Layers layer,
                       const geom::FieldCell<Vect>& fc_corr) override {
    auto sem = mesh.GetSem("corr");
    for (size_t n = 0; n < dim; ++n) {
      if (sem.Nested("dir-corr")) {
        v_solver_[n]->CorrectField(layer, GetComponent(fc_corr, n));
      }
    }
    if (sem("tovect")) {
      CopyToVector(layer);
    }
  }
  const geom::FieldCell<Expr>& GetVelocityEquations(size_t comp) override {
    return v_solver_[comp]->GetEquations();
  }
  geom::MapFace<std::shared_ptr<ConditionFace>>&
  GetVelocityCond(size_t comp) {
    return v_mf_velocity_cond_[comp];
  }
};

template <class Mesh>
class FluidSolver : public UnsteadyIterativeSolver {
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  geom::FieldCell<Scal>* p_fc_density_;
  geom::FieldCell<Scal>* p_fc_viscosity_;
  geom::FieldCell<Vect>* p_fc_force_;
  geom::FieldCell<Vect>* p_fc_stforce_;
  geom::FieldFace<Vect>* p_ff_stforce_;
  geom::FieldCell<Scal>* p_fc_volume_source_;
  geom::FieldCell<Scal>* p_fc_mass_source_;
  Vect meshvel_;

 public:
  FluidSolver(double time, double time_step,
                      geom::FieldCell<Scal>* p_fc_density,
                      geom::FieldCell<Scal>* p_fc_viscosity,
                      geom::FieldCell<Vect>* p_fc_force,
                      geom::FieldCell<Vect>* p_fc_stforce,
                      geom::FieldFace<Vect>* p_ff_stforce,
                      geom::FieldCell<Scal>* p_fc_volume_source,
                      geom::FieldCell<Scal>* p_fc_mass_source,
                      double convergence_tolerance,
                      size_t num_iterations_limit,
                      Vect meshvel)
      : UnsteadyIterativeSolver(time, time_step,
                                convergence_tolerance, num_iterations_limit)
      , p_fc_density_(p_fc_density)
      , p_fc_viscosity_(p_fc_viscosity)
      , p_fc_force_(p_fc_force)
      , p_fc_stforce_(p_fc_stforce)
      , p_ff_stforce_(p_ff_stforce)
      , p_fc_volume_source_(p_fc_volume_source)
      , p_fc_mass_source_(p_fc_mass_source)
      , meshvel_(meshvel)
  {

  }
  virtual const geom::FieldCell<Vect>& GetVelocity() = 0;
  virtual const geom::FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual const geom::FieldCell<Scal>& GetPressure() = 0;
  virtual const geom::FieldCell<Scal>& GetPressure(Layers layer) = 0;
  virtual const geom::FieldFace<Scal>& GetVolumeFlux() = 0;
  virtual const geom::FieldFace<Scal>& GetVolumeFlux(Layers layer) = 0;
  virtual Vect GetMeshVel() { return meshvel_; }
  virtual void SetMeshVel(Vect meshvel) { meshvel_ = meshvel; }
  virtual double GetAutoTimeStep() { return GetTimeStep(); }
};

class ConditionFaceFluid : public ConditionFace {};

class ConditionCellFluid : public ConditionCell {};

namespace fluid_condition {

template <class Mesh>
class NoSlipWall : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  virtual Vect GetVelocity() const = 0;
};

template <class Mesh>
class NoSlipWallFixed : public NoSlipWall<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;
 public:
  NoSlipWallFixed(Vect velocity)
      : velocity_(velocity)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
};

template <class Mesh>
class Inlet : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  virtual Vect GetVelocity() const = 0;
  virtual size_t GetNeighbourCellId() const = 0;
};

template <class Mesh>
class InletFixed : public Inlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;
  size_t neighbour_cell_id_;

 public:
  InletFixed(Vect velocity, size_t neighbour_cell_id)
      : velocity_(velocity)
      , neighbour_cell_id_(neighbour_cell_id)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  size_t GetNeighbourCellId() const override {
    return neighbour_cell_id_;
  }
};

template <class Mesh>
class Outlet : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect velocity) = 0;
  virtual size_t GetNeighbourCellId() const = 0;
};

template <class Mesh>
class OutletAuto : public Outlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;
  size_t neighbour_cell_id_;

 public:
  OutletAuto(size_t neighbour_cell_id)
      : velocity_(Vect::kZero)
      , neighbour_cell_id_(neighbour_cell_id)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  void SetVelocity(Vect velocity) override {
    velocity_ = velocity;
  }
  size_t GetNeighbourCellId() const override {
    return neighbour_cell_id_;
  }
};

template <class Mesh>
class GivenVelocityAndPressure : public ConditionCellFluid {
  using Vect = typename Mesh::Vect;
  using Scal = typename Mesh::Scal;
 public:
  virtual Vect GetVelocity() const = 0;
  virtual Scal GetPressure() const = 0;
};

template <class Mesh>
class GivenVelocityAndPressureFixed : public GivenVelocityAndPressure<Mesh> {
  using Vect = typename Mesh::Vect;
  using Scal = typename Mesh::Scal;
  Vect velocity_;
  Scal pressure_;

 public:
  GivenVelocityAndPressureFixed(Vect velocity, Scal pressure)
      : velocity_(velocity)
      , pressure_(pressure)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  Scal GetPressure() const override {
    return pressure_;
  }
};

template <class Mesh>
class GivenPressure : public ConditionCellFluid {
  using Scal = typename Mesh::Scal;
 public:
  virtual Scal GetPressure() const = 0;
};

template <class Mesh>
class GivenPressureFixed : public GivenPressure<Mesh> {
  using Scal = typename Mesh::Scal;
  Scal pressure_;

 public:
  GivenPressureFixed(Scal pressure)
      : pressure_(pressure)
  {}
  Scal GetPressure() const override {
    return pressure_;
  }
};

} // namespace fluid_condition


template <class Mesh>
std::shared_ptr<ConditionFaceFluid> Parse(std::string argstr,
                                          geom::IdxFace idxface,
                                          const Mesh& mesh) {
  using namespace fluid_condition;
  std::stringstream arg(argstr);

  std::string name;
  arg >> name;

  if (name == "wall") {
    typename Mesh::Vect vel;
    arg >> vel;
    return std::make_shared<NoSlipWallFixed<Mesh>>(vel);
  } else if (name == "inlet") {
    typename Mesh::Vect vel;
    arg >> vel;
    return std::make_shared<InletFixed<Mesh>>(
        vel, mesh.GetValidNeighbourCellId(idxface));
  } else if (name == "outlet") {
    return std::make_shared<OutletAuto<Mesh>>(
        mesh.GetValidNeighbourCellId(idxface));
  } else {
    throw std::runtime_error("Parse: Unknown boundary condition type");
  }
}


// TODO: Second order time step
// TODO: Extrapolation for first iteration
template <class Mesh>
class FluidSimple : public FluidSolver<Mesh> {
  Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  geom::FieldCell<Vect> fc_force_;
  Scal velocity_relaxation_factor_;
  Scal pressure_relaxation_factor_;
  Scal rhie_chow_factor_;
  LayersData<geom::FieldFace<Scal>> ff_vol_flux_;
  std::shared_ptr<ConvectionDiffusionImplicit<Mesh>> conv_diff_solver_;

  LayersData<geom::FieldCell<Scal>> fc_pressure_;
  geom::FieldCell<Scal> fc_kinematic_viscosity_;
  geom::FieldFace<Scal> ff_kinematic_viscosity_;

  geom::MapFace<std::shared_ptr<ConditionFaceFluid>> mf_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_velocity_cond_;
  // TODO: Const specifier for ConditionFace*

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_pressure_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_pressure_grad_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_force_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_pressure_corr_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_viscosity_cond_;

  geom::MapCell<std::shared_ptr<ConditionCellFluid>> mc_cond_;
  geom::MapCell<std::shared_ptr<ConditionCell>> mc_pressure_cond_;
  geom::MapCell<std::shared_ptr<ConditionCell>> mc_velocity_cond_;

  std::shared_ptr<LinearSolver<Scal, IdxCell, Expr>> linear_;

  geom::FieldFace<bool> is_boundary_;

  // common buffers
  geom::FieldFace<Scal> ff_pressure_;
  geom::FieldCell<Vect> fc_pressure_grad_;
  geom::FieldFace<Vect> ff_pressure_grad_;
  geom::FieldCell<Vect> fc_velocity_asterisk_;
  geom::FieldFace<Vect> ff_velocity_asterisk_;
  geom::FieldFace<Scal> ff_volume_flux_asterisk_;
  geom::FieldFace<Scal> ff_volume_flux_interpolated_;
  geom::FieldCell<Scal> fc_diag_coeff_;
  geom::FieldFace<Scal> ff_diag_coeff_;
  geom::FieldFace<Expr> ff_volume_flux_corr_;
  geom::FieldCell<Expr> fc_pressure_corr_system_;
  geom::FieldCell<Scal> fc_pressure_corr_;
  geom::FieldCell<Vect> fc_pressure_corr_grad_;
  geom::FieldFace<Vect> ff_ext_force_;
  geom::FieldCell<Vect> fc_ext_force_restored_;
  geom::FieldFace<Vect> ff_ext_force_restored_;
  geom::FieldFace<Vect> ff_stforce_restored_;
  geom::FieldCell<Vect> fc_velocity_corr_;
  geom::FieldCell<Scal> fc_velcomp_;
  // / needed for MMIM, now disabled
  //geom::FieldFace<Vect> ff_velocity_iter_prev_;
  // LS
  std::vector<Scal> lsa_, lsb_, lsx_;

  MultiTimer<std::string>* timer_;
  bool time_second_order_;
  bool simpler_;
  bool force_geometric_average_;
  Scal guess_extrapolation_;

  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    for (auto it = mf_cond_.cbegin();
        it != mf_cond_.cend(); ++it) {
      IdxFace idxface = it->GetIdx();
      ConditionFaceFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<NoSlipWall<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionFaceValueFixed<Vect>*>(
            mf_velocity_cond_[idxface].get()) =
            ConditionFaceValueFixed<Vect>(cond->GetVelocity());
      } else if (auto cond = dynamic_cast<Inlet<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionFaceValueFixed<Vect>*>(
            mf_velocity_cond_[idxface].get()) =
            ConditionFaceValueFixed<Vect>(cond->GetVelocity());
      } else if (auto cond = dynamic_cast<Outlet<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionFaceValueFixed<Vect>*>(
            mf_velocity_cond_[idxface].get()) =
            ConditionFaceValueFixed<Vect>(cond->GetVelocity());
      } else {
        throw std::runtime_error("Unknown fluid condition");
      }
    }


    for (auto it = mc_cond_.cbegin();
        it != mc_cond_.cend(); ++it) {
      IdxCell idxcell = it->GetIdx();
      ConditionCellFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<GivenPressure<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mc_pressure_cond_[idxcell].get()) =
                ConditionCellValueFixed<Scal>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionCellValueFixed<Vect>*>(
            mc_velocity_cond_[idxcell].get()) =
                ConditionCellValueFixed<Vect>(cond->GetVelocity());
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mc_pressure_cond_[idxcell].get()) =
                ConditionCellValueFixed<Scal>(cond->GetPressure());
      } else {
        throw std::runtime_error("Unknown fluid cell condition");
      }
    }
  }
  // TODO: Consider seperate channels in one domain
  void UpdateOutletBaseConditions() {
    auto& m = mesh;
    using namespace fluid_condition;
    // Extrapolate velocity on outlet faces from cell centers
    // and calculate total inlet and outlet volumetric fluxes
    Scal fi = 0.;   // Both should be positive
    Scal fo = 0.;
    Scal ao = 0.;
    for (auto it = mf_cond_.cbegin();
        it != mf_cond_.cend(); ++it) {
      IdxFace i = it->GetIdx();
      ConditionFaceFluid* cb = it->GetValue().get(); // cond base

      if (auto cd = dynamic_cast<Outlet<Mesh>*>(cb)) {
        size_t id = cd->GetNeighbourCellId();
        IdxCell c = m.GetNeighbourCell(i, id);
        Scal w = (id == 0 ? 1. : -1.);
        cd->SetVelocity(this->GetVelocity(Layers::iter_curr)[c]);
        fo += cd->GetVelocity().dot(m.GetSurface(i)) * w;
        ao += m.GetArea(i);
      } else if (auto cd = dynamic_cast<Inlet<Mesh>*>(cb)) {
        size_t id = cd->GetNeighbourCellId();
        Scal w = (id == 0 ? -1. : 1.);
        fi += cd->GetVelocity().dot(m.GetSurface(i)) * w;
      }
    }

    for (auto i : m.Cells()) {
      fi += (*this->p_fc_volume_source_)[i] * m.GetVolume(i);
    }

    Scal velcor = (fi - fo) / ao; // Additive correction for velocity

    // Apply correction on outlet faces
    for (auto it = mf_cond_.cbegin();
        it != mf_cond_.cend(); ++it) {
      IdxFace i = it->GetIdx();
      ConditionFaceFluid* cb = it->GetValue().get(); // cond base

      if (auto cd = dynamic_cast<Outlet<Mesh>*>(cb)) {
        size_t id = cd->GetNeighbourCellId();
        Scal w = (id == 0 ? 1. : -1.);
        Vect n = m.GetNormal(i);
        cd->SetVelocity(cd->GetVelocity() + n * (velcor * w));
      }
    }

    std::cerr 
        << "fi = " << fi 
        << ", fo = " << fo
        << ", ao = " << ao
        << std::endl;
  }

  void CalcExtForce() {
    timer_->Push("fluid.1.force-correction");
    // Interpolate force to faces (considered as given force)
    ff_ext_force_ = Interpolate(
        *this->p_fc_force_, mf_force_cond_, mesh,
        force_geometric_average_);
    ff_stforce_restored_ = Interpolate(
        *this->p_fc_stforce_, mf_force_cond_, mesh,
        force_geometric_average_);
    fc_ext_force_restored_.Reinit(mesh);

    for (auto idxcell : mesh.SuCells()) {
      Vect sum = Vect::kZero;
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
        sum += mesh.GetSurface(idxface) *
            ff_ext_force_[idxface].dot(mesh.GetNormal(idxface)) *
            mesh.GetCenter(idxcell).dist(mesh.GetCenter(idxface));
      }
      fc_ext_force_restored_[idxcell] = sum / mesh.GetVolume(idxcell);
    }
    // Interpolated restored force to faces (needed later)
    ff_ext_force_restored_ = Interpolate(
        fc_ext_force_restored_, mf_force_cond_, mesh,
        force_geometric_average_);
    timer_->Pop();
  }
  void CalcKinematicViscosity() {
    fc_kinematic_viscosity_.Reinit(mesh);
    for (auto idxcell : mesh.AllCells()) {
      fc_kinematic_viscosity_[idxcell] =
          (*this->p_fc_viscosity_)[idxcell];
    }
    ff_kinematic_viscosity_ = Interpolate(
        fc_kinematic_viscosity_, mf_viscosity_cond_, mesh,
        force_geometric_average_);
  }

 public:
  // TODO: Add Comm for initial fields or require taht from user.
  FluidSimple(Mesh& mesh,
              const geom::FieldCell<Vect>& fc_velocity_initial,
              const geom::MapFace<std::shared_ptr<ConditionFaceFluid>>&
              mf_cond,
              const geom::MapCell<std::shared_ptr<ConditionCellFluid>>&
              mc_cond,
              Scal velocity_relaxation_factor,
              Scal pressure_relaxation_factor,
              Scal rhie_chow_factor,
              geom::FieldCell<Scal>* p_fc_density,
              geom::FieldCell<Scal>* p_fc_viscosity,
              geom::FieldCell<Vect>* p_fc_force,
              geom::FieldCell<Vect>* p_fc_stforce,
              geom::FieldFace<Vect>* p_ff_stforce,
              geom::FieldCell<Scal>* p_fc_volume_source,
              geom::FieldCell<Scal>* p_fc_mass_source,
              double time, double time_step,
              const LinearSolverFactory& linear_factory_velocity,
              const LinearSolverFactory& linear_factory_pressure,
              double convergence_tolerance,
              size_t num_iterations_limit,
              MultiTimer<std::string>* timer,
              bool time_second_order,
              bool simpler,
              bool force_geometric_average,
              Scal guess_extrapolation = 0.,
              Vect meshvel=0)
      : FluidSolver<Mesh>(time, time_step, p_fc_density, p_fc_viscosity,
                    p_fc_force, p_fc_stforce, p_ff_stforce,
                    p_fc_volume_source, p_fc_mass_source,
                    convergence_tolerance, num_iterations_limit, meshvel)
      , mesh(mesh)
      , fc_force_(mesh)
      , velocity_relaxation_factor_(velocity_relaxation_factor)
      , pressure_relaxation_factor_(pressure_relaxation_factor)
      , rhie_chow_factor_(rhie_chow_factor)
      , mf_cond_(mf_cond)
      , mc_cond_(mc_cond)
      , ff_volume_flux_corr_(mesh)
      , fc_pressure_corr_system_(mesh)
      , timer_(timer)
      , time_second_order_(time_second_order)
      , simpler_(simpler)
      , force_geometric_average_(force_geometric_average)
      , guess_extrapolation_(guess_extrapolation)
  {
    linear_ = linear_factory_pressure.Create<Scal, IdxCell, Expr>();

    using namespace fluid_condition;

    is_boundary_.Reinit(mesh, false);
    for (auto it = mf_cond_.cbegin(); it != mf_cond_.cend(); ++it) {
      IdxFace idxface = it->GetIdx();
      is_boundary_[idxface] = true;
      ConditionFaceFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<NoSlipWall<Mesh>*>(cond_generic)) {
        mf_velocity_cond_[idxface] =
            std::make_shared<
            ConditionFaceValueFixed<Vect>>(cond->GetVelocity());
        mf_pressure_cond_[idxface] =
            std::make_shared<ConditionFaceExtrapolation>();
      } else if (auto cond = dynamic_cast<Inlet<Mesh>*>(cond_generic)) {
        mf_velocity_cond_[idxface] =
            std::make_shared<
            ConditionFaceValueFixed<Vect>>(cond->GetVelocity());
        mf_pressure_cond_[idxface] =
            std::make_shared<ConditionFaceExtrapolation>();
      } else if (auto cond = dynamic_cast<Outlet<Mesh>*>(cond_generic)) {
        mf_velocity_cond_[idxface] =
            std::make_shared<
            ConditionFaceValueFixed<Vect>>(cond->GetVelocity());
        mf_pressure_cond_[idxface] =
            std::make_shared<ConditionFaceExtrapolation>();
      } else {
        throw std::runtime_error("Unknown fluid condition");
      }

      mf_pressure_grad_cond_[idxface] =
          std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect::kZero);
      mf_force_cond_[idxface] =
          std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect::kZero);
      mf_pressure_corr_cond_[idxface] =
          std::make_shared<ConditionFaceExtrapolation>();
      mf_viscosity_cond_[idxface] =
          std::make_shared<ConditionFaceDerivativeFixed<Scal>>(0.);
    }

    for (auto it = mc_cond_.cbegin();
        it != mc_cond_.cend(); ++it) {
      IdxCell idxcell = it->GetIdx();
      ConditionCellFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<GivenPressure<Mesh>*>(cond_generic)) {
        mc_pressure_cond_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Scal>>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<Mesh>*>(cond_generic)) {
        mc_pressure_cond_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Scal>>(cond->GetPressure());
        mc_velocity_cond_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Vect>>(cond->GetVelocity());
      } else {
        throw std::runtime_error("Unknown fluid cell condition");
      }
    }

    conv_diff_solver_ = std::make_shared<
        ConvectionDiffusionImplicit<Mesh>>(
            mesh, fc_velocity_initial,
            mf_velocity_cond_, mc_velocity_cond_,
            velocity_relaxation_factor_,
            p_fc_density, &ff_kinematic_viscosity_, &fc_force_,
            &ff_vol_flux_.iter_prev,
            time, time_step,
            linear_factory_velocity,
            convergence_tolerance, num_iterations_limit,
            time_second_order_, guess_extrapolation_);

    fc_pressure_.time_curr.Reinit(mesh, 0.);
    fc_pressure_.time_prev = fc_pressure_.time_curr;

    // Calc initial volume fluxes
    fc_velocity_asterisk_ = conv_diff_solver_->GetVelocity();
    ff_velocity_asterisk_ = Interpolate(
        fc_velocity_asterisk_, mf_velocity_cond_, mesh);
    ff_vol_flux_.time_curr.Reinit(mesh, 0.);
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.time_curr[idxface] =
          ff_velocity_asterisk_[idxface].dot(mesh.GetSurface(idxface));
    }
    // Apply meshvel
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.time_curr[idxface] -= 
          this->meshvel_.dot(mesh.GetSurface(idxface));
    }

    ff_vol_flux_.time_prev = ff_vol_flux_.time_curr;
    ff_volume_flux_interpolated_.Reinit(mesh, 0.);

    ff_pressure_ = Interpolate(fc_pressure_.time_curr, mf_pressure_cond_, mesh);
    fc_pressure_grad_ = Gradient(ff_pressure_, mesh);
    ff_pressure_grad_ = Interpolate(
        fc_pressure_grad_, mf_pressure_grad_cond_, mesh);
  }
  void StartStep() override {
    auto sem = mesh.GetSem("fluid-start");
    if (sem("convdiff-init")) {
      this->ClearIterationCount();
      if (IsNan(fc_pressure_.time_curr)) {
        throw std::runtime_error("NaN initial pressure");
      }
      conv_diff_solver_->SetTimeStep(this->GetTimeStep());
    }

    if (sem.Nested("convdiff-start")) {
      conv_diff_solver_->StartStep();
    }

    if (sem("convdiff-start")) {
      fc_pressure_.iter_curr = fc_pressure_.time_curr;
      ff_vol_flux_.iter_curr = ff_vol_flux_.time_curr;
      for (auto idxcell : mesh.SuCells()) {
        fc_pressure_.iter_curr[idxcell] +=
            (fc_pressure_.time_curr[idxcell] - fc_pressure_.time_prev[idxcell]) *
            guess_extrapolation_;
      }
      for (auto idxface : mesh.Faces()) {
        ff_vol_flux_.iter_curr[idxface] +=
            (ff_vol_flux_.time_curr[idxface] - ff_vol_flux_.time_prev[idxface]) *
            guess_extrapolation_;
      }
    }
  }
  // TODO: rewrite norm() using dist() where needed
  void MakeIteration() override {
    auto sem = mesh.GetSem("fluid-iter");
    auto& m = mesh;

    auto& fc_pressure_prev = fc_pressure_.iter_prev;
    auto& fc_pressure_curr = fc_pressure_.iter_curr;

    if (sem("pgrad")) {
      fc_pressure_prev = fc_pressure_curr;
      ff_vol_flux_.iter_prev = ff_vol_flux_.iter_curr;

      UpdateOutletBaseConditions();
      UpdateDerivedConditions();

      CalcExtForce();

      CalcKinematicViscosity();

      timer_->Push("fluid.0.pressure-gradient");
      ff_pressure_ = Interpolate(fc_pressure_prev, mf_pressure_cond_, mesh);
      fc_pressure_grad_ = Gradient(ff_pressure_, mesh);
      ff_pressure_grad_ = Interpolate(
          fc_pressure_grad_, mf_pressure_grad_cond_, mesh);
      timer_->Pop();

      // initialize force with zero
      fc_force_.Reinit(mesh, Vect(0));
    }
      // append viscous term
    if (sem("explvisc")) {
      timer_->Push("fluid.1a.explicit-viscosity");
      for (size_t n = 0; n < dim; ++n) {
        fc_velcomp_ = GetComponent(
            conv_diff_solver_->GetVelocity(Layers::iter_curr), n);
        auto ff = Interpolate(fc_velcomp_, 
                              conv_diff_solver_->GetVelocityCond(n), mesh);
        auto gc = Gradient(ff, mesh);
        auto gf = Interpolate(gc, mf_force_cond_, mesh); // adhoc: zero-der cond
        for (auto idxcell : mesh.SuCells()) {
          Vect sum = Vect::kZero;
          for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
            IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
            sum += gf[idxface] * (ff_kinematic_viscosity_[idxface] * 
                mesh.GetOutwardSurface(idxcell, i)[n]);
          }
          fc_force_[idxcell] += sum / mesh.GetVolume(idxcell);
        }
      }
      timer_->Pop();
    }

    if (sem("forceappend")) {
      // append to force
      for (auto idxcell : mesh.AllCells()) {
        fc_force_[idxcell] +=
            fc_pressure_grad_[idxcell] * (-1.) +
            fc_ext_force_restored_[idxcell] +
            (*this->p_fc_stforce_)[idxcell] +
            // Volume source momentum compensation:
            conv_diff_solver_->GetVelocity(Layers::iter_curr)[idxcell] *
            ((*this->p_fc_density_)[idxcell] *
            (*this->p_fc_volume_source_)[idxcell] -
            (*this->p_fc_mass_source_)[idxcell]);
      }

      timer_->Push("fluid.2.convection-diffusion");
    }

    if (sem.Nested("convdiff-iter")) {
      conv_diff_solver_->MakeIteration();
    }

    if (sem("diag-comm")) {
      timer_->Pop();

      fc_diag_coeff_.Reinit(mesh);
      for (auto idxcell : mesh.Cells()) {
        Scal sum = 0.;
        for (size_t n = 0; n < dim; ++n) {
          sum += conv_diff_solver_->GetVelocityEquations(n)[idxcell].CoeffSum();
        }
        fc_diag_coeff_[idxcell] = sum / dim;
      }

      mesh.Comm(&fc_diag_coeff_);
    }
      
    if (sem("pcorr-assemble")) {
      // diag coeff condition
      geom::MapFace<std::shared_ptr<ConditionFace>> dcc;
      for (auto i : mesh.Faces()) {
        if (is_boundary_[i]) {
          dcc[i] = std::make_shared
              <solver::ConditionFaceDerivativeFixed<Scal>>(Scal(0));
        }
      }

      // Define ff_diag_coeff_ on inner faces only
      ff_diag_coeff_ = Interpolate(fc_diag_coeff_, dcc, mesh);

      fc_velocity_asterisk_ = conv_diff_solver_->GetVelocity(Layers::iter_curr);

      ff_velocity_asterisk_ = Interpolate(
          fc_velocity_asterisk_, mf_velocity_cond_, mesh);

      // // needed for MMIM, now disabled
      //ff_velocity_iter_prev_ = Interpolate(
      //    conv_diff_solver_->GetVelocity(Layers::iter_prev),
      //    mf_velocity_cond_, mesh);

      // Calc volumetric flux (asterisk)
      // using momentum interpolation (Rhie-Chow)
      // including hydrostatics-correction
      // TODO: Extend hydrostatics-correction on a non-uniform mesh
      timer_->Push("fluid.3.momentum-interpolation");
      ff_volume_flux_asterisk_.Reinit(mesh);
      for (auto idxface : mesh.Faces()) {
        const auto volume_flux_interpolated =
            ff_velocity_asterisk_[idxface].dot(mesh.GetSurface(idxface));
        if (!is_boundary_[idxface]) {
          IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
          IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
          Vect dm = mesh.GetVectToCell(idxface, 0);
          Vect dp = mesh.GetVectToCell(idxface, 1);
          const auto pressure_surface_derivative_wide =
              (ff_pressure_grad_[idxface] -
              //ff_stforce_restored_[idxface] -
              ff_ext_force_restored_[idxface]).dot(mesh.GetSurface(idxface));
          const auto pressure_surface_derivative_compact =
              (fc_pressure_prev[cp] - fc_pressure_prev[cm]) /
              (dp - dm).norm() * mesh.GetArea(idxface) -
              //(*this->p_ff_stforce_)[idxface].dot(mesh.GetNormal(idxface)) - 
              ff_ext_force_[idxface].dot(mesh.GetSurface(idxface));
          //const auto mmim = (1. - velocity_relaxation_factor_) *
          //    (ff_vol_flux_.iter_prev[idxface] -
          //     ff_velocity_iter_prev_[idxface].dot(mesh.GetSurface(idxface)));
          ff_volume_flux_asterisk_[idxface] =
              volume_flux_interpolated +
              rhie_chow_factor_ * (pressure_surface_derivative_wide -
              pressure_surface_derivative_compact) / ff_diag_coeff_[idxface] +
              0; //mmim; // TODO: Test MMIM
        } else {
          ff_volume_flux_asterisk_[idxface] =
              ff_velocity_asterisk_[idxface].dot(mesh.GetSurface(idxface));
        }
      }

      // Apply meshvel
      for (auto idxface : mesh.Faces()) {
        ff_volume_flux_asterisk_[idxface] -= 
            this->meshvel_.dot(mesh.GetSurface(idxface));
      }

      timer_->Pop();

      // TODO: Rename SurfaceVelocity to MassFlux or VolumeFlux

      timer_->Push("fluid.4.volume-flux");
      // TODO: Rename velocity_corr to smth
      // (it's actually full velocity not just correction)
      // (same for ff_volume_flux_corr_)
      for (auto idxface : mesh.Faces()) {
        auto& expr = ff_volume_flux_corr_[idxface];
        expr.Clear();
        IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
        IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
        Vect dm = mesh.GetVectToCell(idxface, 0);
        Vect dp = mesh.GetVectToCell(idxface, 1);
        auto coeff = - mesh.GetArea(idxface) /
            ((dp - dm).norm() * ff_diag_coeff_[idxface]);
        if (is_boundary_[idxface]) {
          coeff = 0.;
        }
        expr.InsertTerm(-coeff, cm);
        expr.InsertTerm(coeff, cp);
        // adhoc for periodic
        expr.SortTerms(true);
        expr.SetConstant(ff_volume_flux_asterisk_[idxface]);
      }
      timer_->Pop();

      timer_->Push("fluid.5.pressure-system");
      for (auto idxcell : mesh.Cells()) {
        auto& eqn = fc_pressure_corr_system_[idxcell];
        Expr flux_sum;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          flux_sum +=
              ff_volume_flux_corr_[idxface] *
              mesh.GetOutwardFactor(idxcell, i);
        }
        eqn =
            flux_sum -
            Expr((*this->p_fc_volume_source_)[idxcell] *
                 mesh.GetVolume(idxcell));
      }

      // Account for cell conditions for pressure
      for (auto it = mc_pressure_cond_.cbegin();
          it != mc_pressure_cond_.cend(); ++it) {
        IdxCell idxcell(it->GetIdx());
        ConditionCell* cond = it->GetValue().get();
        if (auto cond_value = dynamic_cast<ConditionCellValue<Scal>*>(cond)) {
          for (auto idxlocal : mesh.Cells()) {
            auto& eqn = fc_pressure_corr_system_[idxlocal];
            if (idxlocal == idxcell) { 
              eqn.SetKnownValueDiag(idxcell, cond_value->GetValue());
            } else {
              // Substitute value to obtain symmetrix matrix
              eqn.SetKnownValue(idxcell, cond_value->GetValue());
            }
          }
        }
      }
    }

    if (sem("pcorr-solve")) {
      timer_->Pop();
      auto l = ConvertLs(fc_pressure_corr_system_, lsa_, lsb_, lsx_, mesh);
      using T = typename Mesh::LS::T;
      l.t = T::symm;
      m.Solve(l);
      timer_->Push("fluid.6.pressure-solve");
    }

    if (sem("pcorr-comm")) {
      timer_->Pop();

      timer_->Push("fluid.7.correction");

      // Copy solution
      fc_pressure_corr_.Reinit(mesh);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc_pressure_corr_[c] = lsx_[i++];
      }
      
      // Comm pressure correction (needed for flux correction)
      m.Comm(&fc_pressure_corr_);
    }

    if (sem("pcorr-apply")) {
      // Correct pressure
      for (auto idxcell : mesh.Cells()) {
        fc_pressure_curr[idxcell] = fc_pressure_prev[idxcell] +
            pressure_relaxation_factor_ * fc_pressure_corr_[idxcell];
      }
      m.Comm(&fc_pressure_curr);

      fc_pressure_corr_grad_ = Gradient(
          Interpolate(fc_pressure_corr_, mf_pressure_corr_cond_, mesh),
          mesh);

      // Correct the velocity
      fc_velocity_corr_.Reinit(mesh);
      auto& u = conv_diff_solver_->GetVelocity(Layers::iter_curr);
      for (auto idxcell : mesh.Cells()) {
        fc_velocity_corr_[idxcell] =
        //    u[idxcell] * (-1); // XXX: zero velocity
            fc_pressure_corr_grad_[idxcell] / (-fc_diag_coeff_[idxcell]);
      }
    }

    if (sem.Nested("convdiff-corr")) {
      // correct and comm
      conv_diff_solver_->CorrectVelocity(Layers::iter_curr, fc_velocity_corr_);
    }

    if (sem("pcorr-fluxes")) {
      // Calc divergence-free volume fluxes
      for (auto idxface : mesh.Faces()) {
        ff_vol_flux_.iter_curr[idxface] =
            ff_volume_flux_corr_[idxface].Evaluate(fc_pressure_corr_);
      }
      //ff_vol_flux_.iter_curr.Reinit(mesh, 0); // XXX: zero velocity
      timer_->Pop();

      // TODO: SIMPLER removed

      this->IncIterationCount();
      // m.Comm(&fc_velocity_curr);  // Comm done by CorrectVelocity
    }
  }
  void FinishStep() override {
    auto sem = mesh.GetSem("fluid-finish");
    if (sem("inctime")) {
      fc_pressure_.time_prev = fc_pressure_.time_curr;
      ff_vol_flux_.time_prev = ff_vol_flux_.time_curr;
      fc_pressure_.time_curr = fc_pressure_.iter_curr;
      ff_vol_flux_.time_curr = ff_vol_flux_.iter_curr;
      if (IsNan(fc_pressure_.time_curr)) {
        throw std::runtime_error("NaN pressure");
      }
      this->IncTime();
    }
    if (sem.Nested("convdiff-finish")) {
      conv_diff_solver_->FinishStep();
    }
  }
  double GetConvergenceIndicator() const override {
    return conv_diff_solver_->GetConvergenceIndicator();
  }
  const geom::FieldCell<Vect>& GetVelocity() override {
    return conv_diff_solver_->GetVelocity();
  }
  const geom::FieldCell<Vect>& GetVelocity(Layers layer) override {
    return conv_diff_solver_->GetVelocity(layer);
  }
  const geom::FieldCell<Scal>& GetPressure() override {
    return fc_pressure_.time_curr;
  }
  const geom::FieldCell<Scal>& GetPressure(Layers layer) override {
    return fc_pressure_.Get(layer);
  }
  const geom::FieldFace<Scal>& GetVolumeFlux() override {
    return ff_vol_flux_.time_curr;
  }
  const geom::FieldFace<Scal>& GetVolumeFlux(Layers layer) override {
    return ff_vol_flux_.Get(layer);
  }
  double GetAutoTimeStep() override { 
    double dt = 1e10;
    auto& flux = ff_vol_flux_.time_curr;
    for (auto idxcell : mesh.Cells()) {
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
        if (flux[idxface] != 0.) {
          dt = std::min<Scal>(
              dt, std::abs(mesh.GetVolume(idxcell) / flux[idxface]));
        }
      }
    }
    return dt; 
  }
};


} // namespace solver

