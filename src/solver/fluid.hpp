#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "solver.hpp"
#include "metrics.hpp"
#include "conv_diff.hpp"

namespace solver {

template <class Mesh>
class ConvectionDiffusion : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  FieldCell<Scal>* p_fc_density_;
  // TODO: rename to dynamic viscosity
  FieldFace<Scal>* p_ff_kinematic_viscosity_; 
  FieldCell<Vect>* p_fc_force_;
  FieldFace<Scal>* p_ff_vol_flux_;

 public:
  ConvectionDiffusion(double time, double time_step,
                      FieldCell<Scal>* p_fc_density,
                      FieldFace<Scal>* p_ff_kinematic_viscosity,
                      FieldCell<Vect>* p_fc_force,
                      FieldFace<Scal>* p_ff_vol_flux
                      )
      : UnsteadyIterativeSolver(time, time_step)
      , p_fc_density_(p_fc_density)
      , p_ff_kinematic_viscosity_(p_ff_kinematic_viscosity)
      , p_fc_force_(p_fc_force)
      , p_ff_vol_flux_(p_ff_vol_flux)
  {

  }
  virtual const FieldCell<Vect>& GetVelocity() = 0;
  virtual const FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual void CorrectVelocity(Layers layer,
                               const FieldCell<Vect>& fc_corr) = 0;
  virtual const FieldCell<Expr>& GetVelocityEquations(size_t comp) = 0;
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
  LayersData<FieldCell<Vect>> fc_velocity_;

  MapFace<std::shared_ptr<ConditionFace>> mf_velocity_cond_;
  MapCell<std::shared_ptr<ConditionCell>> mc_velocity_cond_;
  FieldCell<Vect>* p_fc_force_;

  VectGeneric<MapFace<std::shared_ptr<ConditionFace>>>
  v_mf_velocity_cond_;
  // TODO: Extract scalar CellCondition
  VectGeneric<std::shared_ptr<Solver>>
  v_solver_;
  VectGeneric<FieldCell<Scal>>
  v_fc_force_;

 public:
  using Par = typename Solver::Par;
  Par* par;
  void CopyToVector(Layers layer) {
    fc_velocity_.Get(layer).Reinit(mesh);
    for (size_t n = 0; n < dim; ++n) {
      SetComponent(fc_velocity_.Get(layer), n, v_solver_[n]->GetField(layer));
    }
  }
  ConvectionDiffusionImplicit(
      Mesh& mesh,
      const FieldCell<Vect>& fc_velocity_initial,
      const MapFace<std::shared_ptr<ConditionFace>>&
      mf_velocity_cond,
      const MapCell<std::shared_ptr<ConditionCell>>&
      mc_velocity_cond,
      FieldCell<Scal>* p_fc_density,
      FieldFace<Scal>* p_ff_kinematic_viscosity,
      FieldCell<Vect>* p_fc_force,
      FieldFace<Scal>* p_ff_vol_flux,
      double t, double dt, Par* par)
      : ConvectionDiffusion<Mesh>(t, dt, p_fc_density, 
                                  p_ff_kinematic_viscosity, 
                                  p_fc_force, p_ff_vol_flux)
      , mesh(mesh)
      , mf_velocity_cond_(mf_velocity_cond)
      , mc_velocity_cond_(mc_velocity_cond)
      , p_fc_force_(p_fc_force)
      , par(par)
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
          mesh, 
          GetComponent(fc_velocity_initial, n),
          v_mf_velocity_cond_[n],
          MapCell<std::shared_ptr<ConditionCell>>() /*TODO empty*/,
          p_fc_density, p_ff_kinematic_viscosity,
          &(v_fc_force_[n]), p_ff_vol_flux, t, dt, par);
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
      this->ClearIter();
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
      this->IncIter();
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
  double GetError() const override {
    if (this->GetIter() == 0) {
      return 1.;
    }
    return CalcDiff(fc_velocity_.iter_curr, fc_velocity_.iter_prev, mesh);
  }
  const FieldCell<Vect>& GetVelocity() override {
    return fc_velocity_.time_curr;
  }
  const FieldCell<Vect>& GetVelocity(Layers layer) override {
    return fc_velocity_.Get(layer);
  }
  void CorrectVelocity(Layers layer,
                       const FieldCell<Vect>& fc_corr) override {
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
  const FieldCell<Expr>& GetVelocityEquations(size_t comp) override {
    return v_solver_[comp]->GetEquations();
  }
  MapFace<std::shared_ptr<ConditionFace>>&
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
  FieldCell<Scal>* p_fc_density_;
  FieldCell<Scal>* p_fc_viscosity_;
  FieldCell<Vect>* p_fc_force_;
  FieldCell<Vect>* p_fc_stforce_;
  FieldFace<Vect>* p_ff_stforce_;
  FieldCell<Scal>* p_fc_volume_source_;
  FieldCell<Scal>* p_fc_mass_source_;

 public:
  FluidSolver(double time, double time_step,
              FieldCell<Scal>* p_fc_density,
              FieldCell<Scal>* p_fc_viscosity,
              FieldCell<Vect>* p_fc_force,
              FieldCell<Vect>* p_fc_stforce,
              FieldFace<Vect>* p_ff_stforce,
              FieldCell<Scal>* p_fc_volume_source,
              FieldCell<Scal>* p_fc_mass_source)
      : UnsteadyIterativeSolver(time, time_step)
      , p_fc_density_(p_fc_density)
      , p_fc_viscosity_(p_fc_viscosity)
      , p_fc_force_(p_fc_force)
      , p_fc_stforce_(p_fc_stforce)
      , p_ff_stforce_(p_ff_stforce)
      , p_fc_volume_source_(p_fc_volume_source)
      , p_fc_mass_source_(p_fc_mass_source)
  {

  }
  virtual const FieldCell<Vect>& GetVelocity() = 0;
  virtual const FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual const FieldCell<Scal>& GetPressure() = 0;
  virtual const FieldCell<Scal>& GetPressure(Layers layer) = 0;
  virtual const FieldFace<Scal>& GetVolumeFlux() = 0;
  virtual const FieldFace<Scal>& GetVolumeFlux(Layers layer) = 0;
  virtual double GetAutoTimeStep() { return GetTimeStep(); }
};

class ConditionFaceFluid : public ConditionFace {
 public:
  ConditionFaceFluid(size_t nci) : ConditionFace(nci) {}
};

class ConditionCellFluid : public ConditionCell {};

namespace fluid_condition {

template <class Mesh>
class NoSlipWall : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  NoSlipWall(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
};

template <class Mesh>
class NoSlipWallFixed : public NoSlipWall<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;
 public:
  NoSlipWallFixed(Vect velocity, size_t nci)
      : NoSlipWall<Mesh>(nci), velocity_(velocity)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
};

template <class Mesh>
class Inlet : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  Inlet(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect velocity) = 0;
};

template <class Mesh>
class InletFixed : public Inlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;

 public:
  InletFixed(Vect velocity, size_t nci)
      : Inlet<Mesh>(nci)
      , velocity_(velocity)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  void SetVelocity(Vect velocity) override {
    velocity_ = velocity;
  }
};

template <class Mesh>
class InletFlux : public Inlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect vel_;
  size_t id_;

 public:
  InletFlux(Vect vel, size_t id, size_t nci)
      : Inlet<Mesh>(nci)
      , vel_(vel)
      , id_(id)
  {}
  Vect GetVelocity() const override {
    return vel_;
  }
  void SetVelocity(Vect vel) override {
    vel_ = vel;
  }
  size_t GetId() {
    return id_;
  }
};

template <class Mesh>
class Outlet : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  Outlet(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect velocity) = 0;
};

template <class Mesh>
class OutletAuto : public Outlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;

 public:
  OutletAuto(size_t nci)
      : Outlet<Mesh>(nci)
      , velocity_(Vect::kZero)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  void SetVelocity(Vect velocity) override {
    velocity_ = velocity;
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
                                          IdxFace idxface,
                                          size_t nc, // neighbour cell id
                                          const Mesh& mesh) {
  using namespace fluid_condition;
  using Vect=  typename Mesh::Vect;
  std::stringstream arg(argstr);

  std::string name;
  arg >> name;

  if (name == "wall") {
    // No-slip wall.
    // wall <velocity>
    Vect vel;
    arg >> vel;
    return std::make_shared<NoSlipWallFixed<Mesh>>(vel, nc);
  } else if (name == "inlet") {
    // Fixed velocity inlet.
    // inlet <velocity>
    Vect vel;
    arg >> vel;
    return std::make_shared<InletFixed<Mesh>>(vel, nc);
  } else if (name == "inletflux") {
    // Fixed flux inlet. Flux defined by given velocity is redistributed
    // over all faces with same id.
    // inletflux <velocity> <id>
    Vect vel;
    int id;
    arg >> vel >> id;
    return std::make_shared<InletFlux<Mesh>>(vel, id, nc);
  } else if (name == "outlet") {
    // Outlet. Velocity is extrapolated from neighbour cells and corrected
    // to yield zero total flux over outlet and inlet faces.
    return std::make_shared<OutletAuto<Mesh>>(nc);
  } else {
    throw std::runtime_error("Parse: Unknown boundary condition type");
  }
}


template <class Mesh>
class FluidSimple : public FluidSolver<Mesh> {
  Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = IdxCell;
  using IdxFace = IdxFace;
  using IdxNode = IdxNode;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  FieldCell<Vect> fc_force_;
  LayersData<FieldFace<Scal>> ff_vol_flux_;
  using CD = ConvectionDiffusionImplicit<Mesh>;
  std::shared_ptr<CD> conv_diff_solver_;
  using CDPar = typename CD::Par;

  LayersData<FieldCell<Scal>> fc_pressure_;
  FieldCell<Scal> fc_kinematic_viscosity_;
  FieldFace<Scal> ff_kinematic_viscosity_;

  MapFace<std::shared_ptr<ConditionFaceFluid>> mf_cond_;

  MapFace<std::shared_ptr<ConditionFace>> mf_velocity_cond_;
  // TODO: Const specifier for ConditionFace*

  MapFace<std::shared_ptr<ConditionFace>> mf_pressure_cond_;

  MapFace<std::shared_ptr<ConditionFace>> mf_pressure_grad_cond_;

  MapFace<std::shared_ptr<ConditionFace>> mf_force_cond_;

  // diag coeff condition
  MapFace<std::shared_ptr<ConditionFace>> mf_dcc_;

  MapFace<std::shared_ptr<ConditionFace>> mf_pressure_corr_cond_;

  MapFace<std::shared_ptr<ConditionFace>> mf_viscosity_cond_;

  MapCell<std::shared_ptr<ConditionCellFluid>> mc_cond_;
  MapCell<std::shared_ptr<ConditionCell>> mc_pressure_cond_;
  MapCell<std::shared_ptr<ConditionCell>> mc_velocity_cond_;

  FieldFace<bool> is_boundary_;

  // used by UpdateOutletBaseConditions():
  Scal olfi_; // inlet flux
  Scal olfo_; // outlet flux
  Scal olao_; // outlet area
  // used by UpdateInletFlux():
  std::vector<Scal> ilft_; // target flux
  std::vector<Scal> ilfe_; // extrapolated flux
  std::vector<Scal> ila_; // area

  // common buffers
  FieldFace<Scal> ff_pressure_;
  FieldCell<Vect> fc_pressure_grad_;
  FieldFace<Vect> ff_pressure_grad_;
  FieldCell<Vect> fc_velocity_asterisk_;
  FieldFace<Vect> ff_velocity_asterisk_;
  FieldFace<Scal> ff_volume_flux_asterisk_;
  FieldFace<Scal> ff_volume_flux_interpolated_;
  FieldCell<Scal> fc_diag_coeff_;
  FieldFace<Scal> ff_diag_coeff_;
  FieldFace<Expr> ff_volume_flux_corr_;
  FieldCell<Expr> fc_pressure_corr_system_;
  FieldCell<Scal> fc_pressure_corr_;
  FieldCell<Vect> fc_pressure_corr_grad_;
  FieldFace<Vect> ff_ext_force_;
  FieldCell<Vect> fc_ext_force_restored_;
  FieldFace<Vect> ff_ext_force_restored_;
  FieldFace<Vect> ff_stforce_restored_;
  FieldCell<Vect> fc_velocity_corr_;
  FieldCell<Scal> fc_velcomp_;
  // / needed for MMIM, now disabled
  //FieldFace<Vect> ff_velocity_iter_prev_;
  // LS
  std::vector<Scal> lsa_, lsb_, lsx_;

  MultiTimer<std::string>* timer_;
  CDPar cdpar;

  // TODO: somhow track dependencies to define execution order
  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    // Face conditions
    for (auto it : mf_cond_) {
      IdxFace i = it.GetIdx();
      ConditionFaceFluid* cb = it.GetValue().get();
      auto p = mf_velocity_cond_[i].get();

      if (auto cd = dynamic_cast<NoSlipWall<Mesh>*>(cb)) {
        auto pd = dynamic_cast<ConditionFaceValueFixed<Vect>*>(p);
        pd->Set(cd->GetVelocity());
      } else if (auto cd = dynamic_cast<Inlet<Mesh>*>(cb)) {
        auto pd = dynamic_cast<ConditionFaceValueFixed<Vect>*>(p);
        pd->Set(cd->GetVelocity());
      } else if (auto cd = dynamic_cast<Outlet<Mesh>*>(cb)) {
        auto pd = dynamic_cast<ConditionFaceValueFixed<Vect>*>(p);
        pd->Set(cd->GetVelocity());
      } else {
        throw std::runtime_error("Unknown fluid condition");
      }
    }

    // Cell conditions
    for (auto it : mc_cond_) {
      IdxCell i = it.GetIdx();
      ConditionCellFluid* cb = it.GetValue().get();

      if (auto cond = dynamic_cast<GivenPressure<Mesh>*>(cb)) {
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mc_pressure_cond_[i].get()) =
                ConditionCellValueFixed<Scal>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<Mesh>*>(cb)) {
        *dynamic_cast<ConditionCellValueFixed<Vect>*>(
            mc_velocity_cond_[i].get()) =
                ConditionCellValueFixed<Vect>(cond->GetVelocity());
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mc_pressure_cond_[i].get()) =
                ConditionCellValueFixed<Scal>(cond->GetPressure());
      } else {
        throw std::runtime_error("Unknown fluid cell condition");
      }
    }
  }
  void UpdateInletFlux() {
    auto& m = mesh;
    using namespace fluid_condition;
    size_t& nid = par->inletflux_numid;

    auto sem = m.GetSem("inletflux");
    if (sem("local")) { 
      ilft_.resize(nid);
      ilfe_.resize(nid);
      ila_.resize(nid);

      for (int id = 0; id < nid; ++id) {
        ilft_[id] = 0.;
        ilfe_[id] = 0.;
        ila_[id] = 0.;
      }

      // Extrapolate velocity to inlet from neighbour cells
      // and compute total fluxes
      auto& vel = this->GetVelocity(Layers::iter_curr);
      for (auto it : mf_cond_) {
        IdxFace i = it.GetIdx();
        ConditionFaceFluid* cb = it.GetValue().get(); // cond base

        size_t nci = cb->GetNci();
        IdxCell c = m.GetNeighbourCell(i, nci);
        if (m.IsInner(c)) {
          if (auto cd = dynamic_cast<InletFlux<Mesh>*>(cb)) {
            size_t id = cd->GetId();
            assert(id < ilft_.size());
            Scal w = (nci == 0 ? -1. : 1.);
            // target flux 
            ilft_[id] += cd->GetVelocity().dot(m.GetSurface(i)) * w;
            // extrapolate velocity
            cd->SetVelocity(vel[c]);
            // extrapolated flux
            ilfe_[id] += cd->GetVelocity().dot(m.GetSurface(i)) * w;
            // area
            ila_[id] += m.GetArea(i);
          }
        }
      }
      
      for (int id = 0; id < nid; ++id) {
        m.Reduce(&ilft_[id], "sum");
        m.Reduce(&ilfe_[id], "sum");
        m.Reduce(&ila_[id], "sum");
      }
    }

    if (sem("corr")) {
      for (int id = 0; id < nid; ++id) {
        // Apply additive correction
        Scal dv = (ilft_[id] - ilfe_[id]) / ila_[id];  // velocity
        for (auto it : mf_cond_) {
          IdxFace i = it.GetIdx();
          ConditionFaceFluid* cb = it.GetValue().get(); // cond base

          if (auto cd = dynamic_cast<InletFlux<Mesh>*>(cb)) {
            size_t nci = cd->GetNci();
            Scal w = (nci == 0 ? -1. : 1.);
            Vect n = m.GetNormal(i);
            cd->SetVelocity(cd->GetVelocity() + n * (dv * w));
          }
        }
      }
    }
  }
  // TODO: Consider seperate channels in one domain
  void UpdateOutletBaseConditions() {
    auto& m = mesh;
    using namespace fluid_condition;

    Scal& fi = olfi_; // total inlet volume flux
    Scal& fo = olfo_; // total outlet volume flux
    Scal& ao = olao_; // total outlet area

    auto sem = m.GetSem("outlet");

    if (sem("local")) { 
      fi = 0.;
      fo = 0.;
      ao = 0.;

      // Extrapolate velocity to outlet from neighbour cells,
      // and compute total fluxes
      for (auto it = mf_cond_.cbegin(); it != mf_cond_.cend(); ++it) {
        IdxFace i = it->GetIdx();
        ConditionFaceFluid* cb = it->GetValue().get(); // cond base

        size_t id = cb->GetNci();
        IdxCell c = m.GetNeighbourCell(i, id);
        if (m.IsInner(c)) {
          if (auto cd = dynamic_cast<Outlet<Mesh>*>(cb)) {
            Scal w = (id == 0 ? 1. : -1.);
            cd->SetVelocity(this->GetVelocity(Layers::iter_curr)[c]);
            fo += cd->GetVelocity().dot(m.GetSurface(i)) * w;
            ao += m.GetArea(i);
          } else if (auto cd = dynamic_cast<Inlet<Mesh>*>(cb)) {
            Scal w = (id == 0 ? -1. : 1.);
            fi += cd->GetVelocity().dot(m.GetSurface(i)) * w;
          }
        }
      }
      
      // Append volume source to inlet flux
      for (auto i : m.Cells()) {
        fi += (*this->p_fc_volume_source_)[i] * m.GetVolume(i);
      }

      m.Reduce(&fi, "sum");
      m.Reduce(&fo, "sum");
      m.Reduce(&ao, "sum");
    }

    if (sem("corr")) {
      Scal velcor = (fi - fo) / ao; // Additive correction for velocity

      // Apply correction on outlet faces
      for (auto it = mf_cond_.cbegin(); it != mf_cond_.cend(); ++it) {
        IdxFace i = it->GetIdx();
        ConditionFaceFluid* cb = it->GetValue().get(); // cond base

        if (auto cd = dynamic_cast<Outlet<Mesh>*>(cb)) {
          size_t id = cd->GetNci();
          Scal w = (id == 0 ? 1. : -1.);
          Vect n = m.GetNormal(i);
          cd->SetVelocity(cd->GetVelocity() + n * (velcor * w));
        }
      }
    }
  }

  void CalcExtForce() {
    timer_->Push("fluid.1.force-correction");
    // Interpolate force to faces (considered as given force)
    ff_ext_force_ = Interpolate(
        *this->p_fc_force_, mf_force_cond_, mesh, par->forcegeom);
    ff_stforce_restored_ = Interpolate(
        *this->p_fc_stforce_, mf_force_cond_, mesh, par->forcegeom);
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
        fc_ext_force_restored_, mf_force_cond_, mesh, par->forcegeom);
    timer_->Pop();
  }
  void CalcKinematicViscosity() {
    fc_kinematic_viscosity_.Reinit(mesh);
    for (auto idxcell : mesh.AllCells()) {
      fc_kinematic_viscosity_[idxcell] =
          (*this->p_fc_viscosity_)[idxcell];
    }
    ff_kinematic_viscosity_ = Interpolate(
        fc_kinematic_viscosity_, mf_viscosity_cond_, mesh, par->forcegeom);
  }

 public:
  struct Par {
    Scal vrelax;   // velocity relaxation factor [0,1]
    Scal prelax;   // pressure relaxation factor [0,1]
    Scal rhie;     // Rhie-Chow factor [0,1] (0 disable, 1 full)
    bool second = false; // second order in time
    bool simpler = false; // Use SIMPLER  TODO: implement SIMPLER
    Scal guessextra = 0;  // next iteration extrapolation weight [0,1]
    Vect meshvel = Vect(0);  // relative mesh velocity
    bool forcegeom = false; // geometric average for force
    size_t inletflux_numid = 0; // reduction for id from 0 to numid-1
  };
  Par* par;
  void Update(CDPar& cdpar, const Par& par);
  // TODO: Add Comm for initial fields or require taht from user.
  FluidSimple(Mesh& mesh,
              const FieldCell<Vect>& fc_velocity_initial,
              const MapFace<std::shared_ptr<ConditionFaceFluid>>&
              mf_cond,
              const MapCell<std::shared_ptr<ConditionCellFluid>>&
              mc_cond,
              FieldCell<Scal>* p_fc_density,
              FieldCell<Scal>* p_fc_viscosity,
              FieldCell<Vect>* p_fc_force,
              FieldCell<Vect>* p_fc_stforce,
              FieldFace<Vect>* p_ff_stforce,
              FieldCell<Scal>* p_fc_volume_source,
              FieldCell<Scal>* p_fc_mass_source,
              double time, double time_step,
              MultiTimer<std::string>* timer,
              Par* par
              )
      : FluidSolver<Mesh>(time, time_step, p_fc_density, p_fc_viscosity,
                    p_fc_force, p_fc_stforce, p_ff_stforce,
                    p_fc_volume_source, p_fc_mass_source)
      , mesh(mesh)
      , fc_force_(mesh)
      , mf_cond_(mf_cond)
      , mc_cond_(mc_cond)
      , ff_volume_flux_corr_(mesh)
      , fc_pressure_corr_system_(mesh)
      , timer_(timer)
      , par(par)
  {
    using namespace fluid_condition;

    is_boundary_.Reinit(mesh, false);
    for (auto it : mf_cond_) {
      IdxFace i = it.GetIdx();
      is_boundary_[i] = true;
      ConditionFaceFluid* cb = it.GetValue().get();
      size_t nci = cb->GetNci();

      if (auto cd = dynamic_cast<NoSlipWall<Mesh>*>(cb)) {
        mf_velocity_cond_[i] =
            std::make_shared<ConditionFaceValueFixed<Vect>>(
                cd->GetVelocity(), nci);
        mf_pressure_cond_[i] =
            std::make_shared<ConditionFaceExtrapolation>(nci);
      } else if (auto cd = dynamic_cast<Inlet<Mesh>*>(cb)) {
        mf_velocity_cond_[i] =
            std::make_shared<ConditionFaceValueFixed<Vect>>(
                cd->GetVelocity(), nci);
        mf_pressure_cond_[i] =
            std::make_shared<ConditionFaceExtrapolation>(nci);
      } else if (auto cd = dynamic_cast<Outlet<Mesh>*>(cb)) {
        mf_velocity_cond_[i] =
            std::make_shared<ConditionFaceValueFixed<Vect>>(
                cd->GetVelocity(), nci);
        mf_pressure_cond_[i] = 
            std::make_shared<ConditionFaceExtrapolation>(nci);
      } else {
        throw std::runtime_error("Unknown fluid condition");
      }

      mf_pressure_grad_cond_[i] =
          std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect::kZero, nci);
      mf_force_cond_[i] =
          std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect::kZero, nci);
      mf_pressure_corr_cond_[i] =
          std::make_shared<ConditionFaceExtrapolation>(nci);
      mf_viscosity_cond_[i] =
          std::make_shared<ConditionFaceDerivativeFixed<Scal>>(0., nci);
      mf_dcc_[i] = 
          std::make_shared<ConditionFaceDerivativeFixed<Scal>>(0, nci);
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

    Update(cdpar, *par);

    conv_diff_solver_ = std::make_shared<
        ConvectionDiffusionImplicit<Mesh>>(
            mesh, fc_velocity_initial,
            mf_velocity_cond_, mc_velocity_cond_,
            p_fc_density, &ff_kinematic_viscosity_, 
            &fc_force_, &ff_vol_flux_.iter_prev,
            time, time_step, &cdpar);

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
    const Vect& meshvel = par->meshvel;
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.time_curr[idxface] -= 
          meshvel.dot(mesh.GetSurface(idxface));
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
      this->ClearIter();
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
      const Scal ge = par->guessextra;
      if (ge != 0.) {
        for (auto idxcell : mesh.SuCells()) {
          fc_pressure_.iter_curr[idxcell] +=
              (fc_pressure_.time_curr[idxcell] - 
               fc_pressure_.time_prev[idxcell]) * ge;
        }
        for (auto idxface : mesh.Faces()) {
          ff_vol_flux_.iter_curr[idxface] +=
              (ff_vol_flux_.time_curr[idxface] - 
               ff_vol_flux_.time_prev[idxface]) * ge;
        }
      }
    }
  }
  // TODO: rewrite norm() using dist() where needed
  void MakeIteration() override {
    auto sem = mesh.GetSem("fluid-iter");
    auto& m = mesh;

    auto& fc_pressure_prev = fc_pressure_.iter_prev;
    auto& fc_pressure_curr = fc_pressure_.iter_curr;

    if (sem.Nested("inletflux")) {
      UpdateInletFlux();
    }

    if (sem.Nested("outlet")) {
      UpdateOutletBaseConditions();
    }

    if (sem("pgrad")) {
      Update(cdpar, *par);
      UpdateDerivedConditions();
      fc_pressure_prev = fc_pressure_curr;
      ff_vol_flux_.iter_prev = ff_vol_flux_.iter_curr;

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
        for (auto c : mesh.Cells()) {
          Vect sum = Vect::kZero;
          for (size_t i = 0; i < mesh.GetNumNeighbourFaces(c); ++i) {
            IdxFace f = mesh.GetNeighbourFace(c, i);
            sum += gf[f] * (ff_kinematic_viscosity_[f] * 
                mesh.GetOutwardSurface(c, i)[n]);
          }
          fc_force_[c] += sum / mesh.GetVolume(c);
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
      // Define ff_diag_coeff_ on inner faces only
      ff_diag_coeff_ = Interpolate(fc_diag_coeff_, mf_dcc_, mesh);

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
      const Scal rh = par->rhie;
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
              rh * (pressure_surface_derivative_wide -
              pressure_surface_derivative_compact) / ff_diag_coeff_[idxface] +
              0; //mmim; // TODO: Test MMIM
        } else {
          ff_volume_flux_asterisk_[idxface] =
              ff_velocity_asterisk_[idxface].dot(mesh.GetSurface(idxface));
        }
      }

      // Apply meshvel
      const Vect& meshvel = par->meshvel;
      for (auto idxface : mesh.Faces()) {
        ff_volume_flux_asterisk_[idxface] -= 
            meshvel.dot(mesh.GetSurface(idxface));
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
      Scal pr = par->prelax;
      for (auto idxcell : mesh.Cells()) {
        fc_pressure_curr[idxcell] = fc_pressure_prev[idxcell] +
            pr * fc_pressure_corr_[idxcell];
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

      this->IncIter();
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
  double GetError() const override {
    return conv_diff_solver_->GetError();
  }
  const FieldCell<Vect>& GetVelocity() override {
    return conv_diff_solver_->GetVelocity();
  }
  const FieldCell<Vect>& GetVelocity(Layers layer) override {
    return conv_diff_solver_->GetVelocity(layer);
  }
  const FieldCell<Scal>& GetPressure() override {
    return fc_pressure_.time_curr;
  }
  const FieldCell<Scal>& GetPressure(Layers layer) override {
    return fc_pressure_.Get(layer);
  }
  const FieldFace<Scal>& GetVolumeFlux() override {
    return ff_vol_flux_.time_curr;
  }
  const FieldFace<Scal>& GetVolumeFlux(Layers layer) override {
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

template <class Mesh>
void FluidSimple<Mesh>::Update(CDPar& d, const Par& p) {
  // Update convdiff parameters
  d.relax = p.vrelax;
  d.guessextra = p.guessextra;
  d.second = p.second;
}

} // namespace solver

