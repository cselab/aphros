#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "util/metrics.h"
#include "convdiffvi.h"
#include "fluid.h"

namespace solver {

template <class Mesh>
class FluidSimple : public FluidSolver<Mesh> {
  Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  FieldCell<Vect> fc_vforce_;
  FieldCell<Vect> fc_cdf_;
  LayersData<FieldFace<Scal>> ff_vol_flux_;
  using CD = ConvectionDiffusionImplicit<Mesh>;
  std::shared_ptr<CD> cd_;

  LayersData<FieldCell<Scal>> fcp_; // pressure

  // TODO: Const specifier for ConditionFace*
  
  // Face conditions
  MapFace<std::shared_ptr<ConditionFaceFluid>> mfc_; // fluid cond
  MapFace<std::shared_ptr<ConditionFace>> mfcw_; // velocity cond
  MapFace<std::shared_ptr<ConditionFace>> mfcp_; // pressure cond
  MapFace<std::shared_ptr<ConditionFace>> mfcgp_; // pressure gradient cond
  MapFace<std::shared_ptr<ConditionFace>> mfcf_; // force cond
  MapFace<std::shared_ptr<ConditionFace>> mfck_; // diag coeff cond
  MapFace<std::shared_ptr<ConditionFace>> mfcpc_; // pressure corr cond
  MapFace<std::shared_ptr<ConditionFace>> mfcd_; // dynamic viscosity cond

  // Cell conditions
  MapCell<std::shared_ptr<ConditionCellFluid>> mcc_; // fluid cell cond
  MapCell<std::shared_ptr<ConditionCell>> mccp_; // pressure cell cond
  MapCell<std::shared_ptr<ConditionCell>> mccw_; // velocity cell cond

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

  // notation:
  // p: pressure
  // gp: pressure gradient
  // w: velocity
  // we: predictor velocity
  // v: volume flux
  // ve: predictor volume flux
  
  // Cell fields:
  FieldCell<Vect> fcgp_; // gradient of pressure 
  FieldCell<Vect> fcwe_; // predictor velocity 
  FieldCell<Scal> fck_; // diag coeff of velocity equation 
  FieldCell<Expr> fcpcs_; // pressure correction linear system
  FieldCell<Scal> fcpc_; // pressure correction
  FieldCell<Vect> fcgpc_; // gradient of pressure correction
  FieldCell<Vect> fffcd_; // force for convdiff
  FieldCell<Vect> fcvc_; // velocity correction
  FieldCell<Scal> fcwo_; // one velocitt component
  FieldCell<Scal> fcdk_; // kinematic viscosity

  // Face fields:
  FieldFace<Scal> ffp_; // pressure
  FieldFace<Vect> ffgp_; // gradient of pressure 
  FieldFace<Vect> ffwe_; // predictor velocity 
  FieldFace<Scal> ffve_; // flux predictor volume 
  FieldFace<Scal> ffv_;  // flux volume 
  FieldFace<Scal> ffk_; // diag coeff of velocity equation 
  FieldFace<Expr> ffvc_; // expression for corrected volume flux
  FieldFace<Vect> fff_;  // restored vector force
  FieldFace<Scal> ffdk_; // kinematic viscosity

  // / needed for MMIM, now disabled
  //FieldFace<Vect> ff_velocity_iter_prev_;
  
  // LS
  std::vector<Scal> lsa_, lsb_, lsx_;

  MultiTimer<std::string>* timer_;

  // TODO: somhow track dependencies to define execution order
  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    // Face conditions
    for (auto it : mfc_) {
      IdxFace i = it.GetIdx();
      ConditionFaceFluid* cb = it.GetValue().get();
      auto p = mfcw_[i].get();

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
    for (auto it : mcc_) {
      IdxCell i = it.GetIdx();
      ConditionCellFluid* cb = it.GetValue().get();

      if (auto cond = dynamic_cast<GivenPressure<Mesh>*>(cb)) {
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mccp_[i].get()) =
                ConditionCellValueFixed<Scal>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<Mesh>*>(cb)) {
        *dynamic_cast<ConditionCellValueFixed<Vect>*>(
            mccw_[i].get()) =
                ConditionCellValueFixed<Vect>(cond->GetVelocity());
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mccp_[i].get()) =
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
      for (auto it : mfc_) {
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
        for (auto it : mfc_) {
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
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
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
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
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
    auto& m = mesh;
    timer_->Push("fluid.1.force-correction");

    // Restore vector force from faces
    fc_vforce_.Reinit(m);
    for (auto c : m.SuCells()) {
      Vect s(0);
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetNeighbourFace(c, q);
        s += m.GetSurface(f) *
            ((*this->p_ff_force_)[f] * m.GetCenter(c).dist(m.GetCenter(f)));
      }
      fc_vforce_[c] = s / m.GetVolume(c);
    }
    // Interpolate vector force to faces
    fff_ = Interpolate(fc_vforce_, mfcf_, m);

    timer_->Pop();
  }
  void CalcKinematicViscosity() {
    fcdk_.Reinit(mesh);
    for (auto idxcell : mesh.AllCells()) {
      fcdk_[idxcell] =
          (*this->p_fc_viscosity_)[idxcell];
    }
    ffdk_ = Interpolate(
        fcdk_, mfcd_, mesh, par->forcegeom);
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
  std::shared_ptr<Par> par;
  Par* GetPar() { return par.get(); }
  void Update(typename CD::Par& cdpar, const Par& par);
  // TODO: Add Comm for initial fields or require taht from user.
  FluidSimple(Mesh& mesh,
              const FieldCell<Vect>& fc_velocity_initial,
              const MapFace<std::shared_ptr<ConditionFaceFluid>>&
              mf_cond,
              const MapCell<std::shared_ptr<ConditionCellFluid>>&
              mc_cond,
              FieldCell<Scal>* p_fc_density,
              FieldCell<Scal>* p_fc_viscosity,
              FieldFace<Scal>* p_ff_force,
              FieldCell<Scal>* p_fc_volume_source,
              FieldCell<Scal>* p_fc_mass_source,
              double time, double time_step,
              MultiTimer<std::string>* timer,
              std::shared_ptr<Par> par
              )
      : FluidSolver<Mesh>(time, time_step, p_fc_density, p_fc_viscosity,
                    p_ff_force, p_fc_volume_source, p_fc_mass_source)
      , mesh(mesh)
      , mfc_(mf_cond)
      , mcc_(mc_cond)
      , ffvc_(mesh)
      , fcpcs_(mesh)
      , timer_(timer)
      , par(par)
  {
    using namespace fluid_condition;

    is_boundary_.Reinit(mesh, false);
    for (auto it : mfc_) {
      IdxFace i = it.GetIdx();
      is_boundary_[i] = true;
      ConditionFaceFluid* cb = it.GetValue().get();
      size_t nci = cb->GetNci();

      if (auto cd = dynamic_cast<NoSlipWall<Mesh>*>(cb)) {
        mfcw_[i] =
            std::make_shared<ConditionFaceValueFixed<Vect>>(
                cd->GetVelocity(), nci);
        mfcp_[i] =
            std::make_shared<ConditionFaceExtrapolation>(nci);
      } else if (auto cd = dynamic_cast<Inlet<Mesh>*>(cb)) {
        mfcw_[i] =
            std::make_shared<ConditionFaceValueFixed<Vect>>(
                cd->GetVelocity(), nci);
        mfcp_[i] =
            std::make_shared<ConditionFaceExtrapolation>(nci);
      } else if (auto cd = dynamic_cast<Outlet<Mesh>*>(cb)) {
        mfcw_[i] =
            std::make_shared<ConditionFaceValueFixed<Vect>>(
                cd->GetVelocity(), nci);
        mfcp_[i] = 
            std::make_shared<ConditionFaceExtrapolation>(nci);
      } else {
        throw std::runtime_error("Unknown fluid condition");
      }

      mfcgp_[i] =
          std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect::kZero, nci);
      mfcf_[i] =
          std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect::kZero, nci);
      mfcpc_[i] =
          std::make_shared<ConditionFaceExtrapolation>(nci);
      mfcd_[i] =
          std::make_shared<ConditionFaceDerivativeFixed<Scal>>(0., nci);
      mfck_[i] = 
          std::make_shared<ConditionFaceDerivativeFixed<Scal>>(0, nci);
    }

    for (auto it = mcc_.cbegin();
        it != mcc_.cend(); ++it) {
      IdxCell idxcell = it->GetIdx();
      ConditionCellFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<GivenPressure<Mesh>*>(cond_generic)) {
        mccp_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Scal>>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<Mesh>*>(cond_generic)) {
        mccp_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Scal>>(cond->GetPressure());
        mccw_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Vect>>(cond->GetVelocity());
      } else {
        throw std::runtime_error("Unknown fluid cell condition");
      }
    }

    // Init convdiff solver
    {
      auto p = std::make_shared<typename CD::Par>();
      Update(*p, *par);

      fc_cdf_.Reinit(mesh, Vect(0));
      cd_ = std::make_shared<
          ConvectionDiffusionImplicit<Mesh>>(
              mesh, fc_velocity_initial,
              mfcw_, mccw_,
              p_fc_density, &ffdk_, 
              &fc_cdf_, &ff_vol_flux_.iter_prev, time, time_step, p);
    }

    fcp_.time_curr.Reinit(mesh, 0.);
    fcp_.time_prev = fcp_.time_curr;

    // Calc initial volume fluxes
    fcwe_ = cd_->GetVelocity();
    ffwe_ = Interpolate(
        fcwe_, mfcw_, mesh);
    ff_vol_flux_.time_curr.Reinit(mesh, 0.);
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.time_curr[idxface] =
          ffwe_[idxface].dot(mesh.GetSurface(idxface));
    }
    // Apply meshvel
    const Vect& meshvel = par->meshvel;
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.time_curr[idxface] -= 
          meshvel.dot(mesh.GetSurface(idxface));
    }

    ff_vol_flux_.time_prev = ff_vol_flux_.time_curr;
    ffv_.Reinit(mesh, 0.);

    ffp_ = Interpolate(fcp_.time_curr, mfcp_, mesh);
    fcgp_ = Gradient(ffp_, mesh);
    ffgp_ = Interpolate(
        fcgp_, mfcgp_, mesh);
  }
  void StartStep() override {
    auto sem = mesh.GetSem("fluid-start");
    if (sem("convdiff-init")) {
      this->ClearIter();
      if (IsNan(fcp_.time_curr)) {
        throw std::runtime_error("NaN initial pressure");
      }
      cd_->SetTimeStep(this->GetTimeStep());
    }

    if (sem.Nested("convdiff-start")) {
      cd_->StartStep();
    }

    if (sem("convdiff-start")) {
      fcp_.iter_curr = fcp_.time_curr;
      ff_vol_flux_.iter_curr = ff_vol_flux_.time_curr;
      const Scal ge = par->guessextra;
      if (ge != 0.) {
        for (auto idxcell : mesh.SuCells()) {
          fcp_.iter_curr[idxcell] +=
              (fcp_.time_curr[idxcell] - 
               fcp_.time_prev[idxcell]) * ge;
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

    auto& fcp_prev = fcp_.iter_prev;
    auto& fcp_curr = fcp_.iter_curr;

    if (sem.Nested("inletflux")) {
      UpdateInletFlux();
    }

    if (sem.Nested("outlet")) {
      UpdateOutletBaseConditions();
    }

    if (sem("pgrad")) {
      Update(*cd_->GetPar(), *par);
      UpdateDerivedConditions();
      fcp_prev = fcp_curr;
      ff_vol_flux_.iter_prev = ff_vol_flux_.iter_curr;

      CalcExtForce();

      CalcKinematicViscosity();

      timer_->Push("fluid.0.pressure-gradient");
      ffp_ = Interpolate(fcp_prev, mfcp_, m);
      fcgp_ = Gradient(ffp_, m);
      ffgp_ = Interpolate(
          fcgp_, mfcgp_, m);
      timer_->Pop();

      // initialize force for convdiff
      fc_cdf_.Reinit(m, Vect(0));
    }

    if (sem("explvisc")) {
      timer_->Push("fluid.1a.explicit-viscosity");
      // append viscous term
      for (size_t n = 0; n < dim; ++n) {
        fcwo_ = GetComponent(
            cd_->GetVelocity(Layers::iter_curr), n);
        auto ff = Interpolate(fcwo_, 
                              cd_->GetVelocityCond(n), m);
        auto gc = Gradient(ff, m);
        auto gf = Interpolate(gc, mfcf_, m); // adhoc: zero-der cond
        for (auto c : m.Cells()) {
          Vect s(0);
          for (auto q : m.Nci(c)) {
            IdxFace f = mesh.GetNeighbourFace(c, q);
            s += gf[f] * 
              (ffdk_[f] * m.GetOutwardSurface(c, q)[n]);
          }
          fc_cdf_[c] += s / mesh.GetVolume(c);
        }
      }
      timer_->Pop();
    }

    if (sem("forceappend")) {
      // append pressure gradient and external force
      for (auto c : m.AllCells()) {
        fc_cdf_[c] += fcgp_[c] * (-1.) + fc_vforce_[c];
      }

      timer_->Push("fluid.2.convection-diffusion");
    }

    if (sem.Nested("convdiff-iter")) {
      cd_->MakeIteration();
    }

    if (sem("diag-comm")) {
      timer_->Pop();

      fck_.Reinit(mesh);
      for (auto idxcell : mesh.Cells()) {
        Scal sum = 0.;
        for (size_t n = 0; n < dim; ++n) {
          sum += cd_->GetVelocityEquations(n)[idxcell].CoeffSum();
        }
        fck_[idxcell] = sum / dim;
      }

      mesh.Comm(&fck_);
    }
      
    if (sem("pcorr-assemble")) {
      // Define ffk_ on inner faces only
      ffk_ = Interpolate(fck_, mfck_, mesh);

      fcwe_ = cd_->GetVelocity(Layers::iter_curr);

      ffwe_ = Interpolate(
          fcwe_, mfcw_, mesh);

      // // needed for MMIM, now disabled
      //ff_velocity_iter_prev_ = Interpolate(
      //    cd_->GetVelocity(Layers::iter_prev),
      //    mfcw_, mesh);

      // Calc volumetric flux (asterisk)
      // using momentum interpolation (Rhie-Chow)
      // including hydrostatics-correction
      // TODO: Extend hydrostatics-correction on a non-uniform mesh
      timer_->Push("fluid.3.momentum-interpolation");
      ffve_.Reinit(m);
      const Scal rh = par->rhie;
      for (auto f : m.Faces()) {
        const auto volume_flux_interpolated =
            ffwe_[f].dot(m.GetSurface(f));
        if (!is_boundary_[f]) {
          IdxCell cm = m.GetNeighbourCell(f, 0);
          IdxCell cp = m.GetNeighbourCell(f, 1);
          Vect dm = m.GetVectToCell(f, 0);
          Vect dp = m.GetVectToCell(f, 1);
          const auto wide =
              (ffgp_[f] - fff_[f]).dot(m.GetSurface(f));
          const auto compact =
              ((fcp_prev[cp] - fcp_prev[cm]) /
              (dp - dm).norm() - (*this->p_ff_force_)[f]) * m.GetArea(f);
          ffve_[f] =
              volume_flux_interpolated +
              rh * (wide - compact) / ffk_[f];
        } else {
          ffve_[f] =
              ffwe_[f].dot(m.GetSurface(f));
        }
      }

      // Apply meshvel
      const Vect& meshvel = par->meshvel;
      for (auto f : m.Faces()) {
        ffve_[f] -= meshvel.dot(m.GetSurface(f));
      }

      timer_->Pop();

      // TODO: Rename SurfaceVelocity to MassFlux or VolumeFlux

      timer_->Push("fluid.4.volume-flux");
      // TODO: Rename velocity_corr to smth
      // (it's actually full velocity not just correction)
      // (same for ffvc_)
      for (auto idxface : mesh.Faces()) {
        auto& expr = ffvc_[idxface];
        expr.Clear();
        IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
        IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
        Vect dm = mesh.GetVectToCell(idxface, 0);
        Vect dp = mesh.GetVectToCell(idxface, 1);
        auto coeff = - mesh.GetArea(idxface) /
            ((dp - dm).norm() * ffk_[idxface]);
        if (is_boundary_[idxface]) {
          coeff = 0.;
        }
        expr.InsertTerm(-coeff, cm);
        expr.InsertTerm(coeff, cp);
        // adhoc for periodic
        expr.SortTerms(true);
        expr.SetConstant(ffve_[idxface]);
      }
      timer_->Pop();

      timer_->Push("fluid.5.pressure-system");
      for (auto idxcell : mesh.Cells()) {
        auto& eqn = fcpcs_[idxcell];
        Expr flux_sum;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          flux_sum +=
              ffvc_[idxface] *
              mesh.GetOutwardFactor(idxcell, i);
        }
        eqn =
            flux_sum -
            Expr((*this->p_fc_volume_source_)[idxcell] *
                 mesh.GetVolume(idxcell));
      }

      // Account for cell conditions for pressure
      for (auto it = mccp_.cbegin();
          it != mccp_.cend(); ++it) {
        IdxCell idxcell(it->GetIdx());
        ConditionCell* cond = it->GetValue().get();
        if (auto cond_value = dynamic_cast<ConditionCellValue<Scal>*>(cond)) {
          for (auto idxlocal : mesh.Cells()) {
            auto& eqn = fcpcs_[idxlocal];
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
      auto l = ConvertLs(fcpcs_, lsa_, lsb_, lsx_, mesh);
      using T = typename Mesh::LS::T;
      l.t = T::symm;
      m.Solve(l);
      timer_->Push("fluid.6.pressure-solve");
    }

    if (sem("pcorr-comm")) {
      timer_->Pop();

      timer_->Push("fluid.7.correction");

      // Copy solution
      fcpc_.Reinit(mesh);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fcpc_[c] = lsx_[i++];
      }
      
      // Comm pressure correction (needed for flux correction)
      m.Comm(&fcpc_);
    }

    if (sem("pcorr-apply")) {
      // Correct pressure
      Scal pr = par->prelax;
      for (auto idxcell : mesh.Cells()) {
        fcp_curr[idxcell] = fcp_prev[idxcell] +
            pr * fcpc_[idxcell];
      }
      m.Comm(&fcp_curr);

      fcgpc_ = Gradient(
          Interpolate(fcpc_, mfcpc_, mesh),
          mesh);

      // Correct the velocity
      fcvc_.Reinit(mesh);
      auto& u = cd_->GetVelocity(Layers::iter_curr);
      for (auto idxcell : mesh.Cells()) {
        fcvc_[idxcell] =
        //    u[idxcell] * (-1); // XXX: zero velocity
            fcgpc_[idxcell] / (-fck_[idxcell]);
      }
    }

    if (sem.Nested("convdiff-corr")) {
      // correct and comm
      cd_->CorrectVelocity(Layers::iter_curr, fcvc_);
    }

    if (sem("pcorr-fluxes")) {
      // Calc divergence-free volume fluxes
      for (auto idxface : mesh.Faces()) {
        ff_vol_flux_.iter_curr[idxface] =
            ffvc_[idxface].Evaluate(fcpc_);
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
      fcp_.time_prev = fcp_.time_curr;
      ff_vol_flux_.time_prev = ff_vol_flux_.time_curr;
      fcp_.time_curr = fcp_.iter_curr;
      ff_vol_flux_.time_curr = ff_vol_flux_.iter_curr;
      if (IsNan(fcp_.time_curr)) {
        throw std::runtime_error("NaN pressure");
      }
      this->IncTime();
    }
    if (sem.Nested("convdiff-finish")) {
      cd_->FinishStep();
    }
  }
  double GetError() const override {
    return cd_->GetError();
  }
  const FieldCell<Vect>& GetVelocity() override {
    return cd_->GetVelocity();
  }
  const FieldCell<Vect>& GetVelocity(Layers layer) override {
    return cd_->GetVelocity(layer);
  }
  const FieldCell<Scal>& GetPressure() override {
    return fcp_.time_curr;
  }
  const FieldCell<Scal>& GetPressure(Layers layer) override {
    return fcp_.Get(layer);
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
void FluidSimple<Mesh>::Update(typename CD::Par& d, const Par& p) {
  // Update convdiff parameters
  d.relax = p.vrelax;
  d.guessextra = p.guessextra;
  d.second = p.second;
}

} // namespace solver

