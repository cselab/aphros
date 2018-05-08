#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "util/metrics.h"
#include "convdiffvi.h"
#include "fluid.h"

namespace solver {

template <class M_>
class FluidSimple : public FluidSolver<M_> {
  using M = M_;
  using P = FluidSolver<M>; // parent
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  using CD = ConvectionDiffusionImplicit<M>; // convdiff solver

  using P::fcr_;
  using P::fcd_;
  using P::fffp_;
  using P::fcsv_;
  using P::fcsm_;

  M& m; // mesh

  LayersData<FieldFace<Scal>> ffv_; // volume flux
  LayersData<FieldCell<Scal>> fcp_; // pressure

  std::shared_ptr<CD> cd_;

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

  FieldFace<bool> ffb_; // is boundary

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
  FieldCell<Vect> fcf_;  // restored vector force
  FieldCell<Vect> fcgp_; // gradient of pressure 
  FieldCell<Vect> fcwe_; // predictor velocity 
  FieldCell<Scal> fck_; // diag coeff of velocity equation 
  FieldCell<Expr> fcpcs_; // pressure correction linear system
  FieldCell<Scal> fcpc_; // pressure correction
  FieldCell<Vect> fcgpc_; // gradient of pressure correction
  FieldCell<Vect> fcvc_; // velocity correction
  FieldCell<Scal> fcwo_; // one velocitt component
  FieldCell<Scal> fcdk_; // kinematic viscosity
  FieldCell<Vect> fc_cdf_; // force for convdiff

  // Face fields:
  FieldFace<Scal> ffp_; // pressure
  FieldFace<Vect> ffgp_; // gradient of pressure 
  FieldFace<Vect> ffwe_; // predictor velocity 
  FieldFace<Scal> ffve_; // flux predictor volume 
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

      if (auto cd = dynamic_cast<NoSlipWall<M>*>(cb)) {
        auto pd = dynamic_cast<ConditionFaceValueFixed<Vect>*>(p);
        pd->Set(cd->GetVelocity());
      } else if (auto cd = dynamic_cast<Inlet<M>*>(cb)) {
        auto pd = dynamic_cast<ConditionFaceValueFixed<Vect>*>(p);
        pd->Set(cd->GetVelocity());
      } else if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
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

      if (auto cond = dynamic_cast<GivenPressure<M>*>(cb)) {
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mccp_[i].get()) =
                ConditionCellValueFixed<Scal>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<M>*>(cb)) {
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
          if (auto cd = dynamic_cast<InletFlux<M>*>(cb)) {
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

          if (auto cd = dynamic_cast<InletFlux<M>*>(cb)) {
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
          if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
            Scal w = (id == 0 ? 1. : -1.);
            cd->SetVelocity(this->GetVelocity(Layers::iter_curr)[c]);
            fo += cd->GetVelocity().dot(m.GetSurface(i)) * w;
            ao += m.GetArea(i);
          } else if (auto cd = dynamic_cast<Inlet<M>*>(cb)) {
            Scal w = (id == 0 ? -1. : 1.);
            fi += cd->GetVelocity().dot(m.GetSurface(i)) * w;
          }
        }
      }
      
      // Append volume source to inlet flux
      for (auto i : m.Cells()) {
        fi += (*fcsv_)[i] * m.GetVolume(i);
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

        if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
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

    // Restore vector force from faces
    fcf_.Reinit(m);
    for (auto c : m.SuCells()) {
      Vect s(0);
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetNeighbourFace(c, q);
        s += m.GetSurface(f) *
            ((*fffp_)[f] * m.GetCenter(c).dist(m.GetCenter(f)));
      }
      fcf_[c] = s / m.GetVolume(c);
    }
    // Interpolate vector force to faces
    fff_ = Interpolate(fcf_, mfcf_, m);

    timer_->Pop();
  }
  void CalcKinematicViscosity() {
    fcdk_.Reinit(m);
    for (auto c : m.AllCells()) {
      fcdk_[c] = (*fcd_)[c];
    }
    ffdk_ = Interpolate(fcdk_, mfcd_, m, par->forcegeom);
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
  FluidSimple(M& m,
              const FieldCell<Vect>& fcw,
              const MapFace<std::shared_ptr<ConditionFaceFluid>>& mfc,
              const MapCell<std::shared_ptr<ConditionCellFluid>>& mcc,
              FieldCell<Scal>* fcr, // density
              FieldCell<Scal>* fcd, // dynamic viscosity
              FieldFace<Scal>* fffp, // force projections on faces
              FieldCell<Scal>* fcsv, // volume source
              FieldCell<Scal>* fcsm, // mass source
              double time, double time_step,
              MultiTimer<std::string>* timer,
              std::shared_ptr<Par> par
              )
      : FluidSolver<M>(time, time_step, fcr, fcd, fffp, fcsv, fcsm)
      , m(m) , mfc_(mfc) , mcc_(mcc) , ffvc_(m) , fcpcs_(m)
      , timer_(timer) , par(par)
  {
    using namespace fluid_condition;

    ffb_.Reinit(m, false);
    for (auto it : mfc_) {
      IdxFace i = it.GetIdx();
      ffb_[i] = true;
      ConditionFaceFluid* cb = it.GetValue().get();
      size_t nci = cb->GetNci();

      if (auto cd = dynamic_cast<NoSlipWall<M>*>(cb)) {
        mfcw_[i] =
            std::make_shared<ConditionFaceValueFixed<Vect>>(
                cd->GetVelocity(), nci);
        mfcp_[i] =
            std::make_shared<ConditionFaceExtrapolation>(nci);
      } else if (auto cd = dynamic_cast<Inlet<M>*>(cb)) {
        mfcw_[i] =
            std::make_shared<ConditionFaceValueFixed<Vect>>(
                cd->GetVelocity(), nci);
        mfcp_[i] =
            std::make_shared<ConditionFaceExtrapolation>(nci);
      } else if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
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
      IdxCell c = it->GetIdx();
      ConditionCellFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<GivenPressure<M>*>(cond_generic)) {
        mccp_[c] =
            std::make_shared<
            ConditionCellValueFixed<Scal>>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<M>*>(cond_generic)) {
        mccp_[c] =
            std::make_shared<
            ConditionCellValueFixed<Scal>>(cond->GetPressure());
        mccw_[c] =
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

      fc_cdf_.Reinit(m, Vect(0));
      cd_ = std::make_shared<
          ConvectionDiffusionImplicit<M>>(
              m, fcw,
              mfcw_, mccw_,
              fcr, &ffdk_, 
              &fc_cdf_, &ffv_.iter_prev, time, time_step, p);
    }

    fcp_.time_curr.Reinit(m, 0.);
    fcp_.time_prev = fcp_.time_curr;

    // Calc initial volume fluxes
    fcwe_ = cd_->GetVelocity();
    ffwe_ = Interpolate(
        fcwe_, mfcw_, m);
    ffv_.time_curr.Reinit(m, 0.);
    for (auto f : m.Faces()) {
      ffv_.time_curr[f] =
          ffwe_[f].dot(m.GetSurface(f));
    }
    // Apply meshvel
    const Vect& meshvel = par->meshvel;
    for (auto f : m.Faces()) {
      ffv_.time_curr[f] -= 
          meshvel.dot(m.GetSurface(f));
    }

    ffv_.time_prev = ffv_.time_curr;

    ffp_ = Interpolate(fcp_.time_curr, mfcp_, m);
    fcgp_ = Gradient(ffp_, m);
    ffgp_ = Interpolate(
        fcgp_, mfcgp_, m);
  }
  void StartStep() override {
    auto sem = m.GetSem("fluid-start");
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
      ffv_.iter_curr = ffv_.time_curr;
      const Scal ge = par->guessextra;
      if (ge != 0.) {
        for (auto c : m.SuCells()) {
          fcp_.iter_curr[c] +=
              (fcp_.time_curr[c] - 
               fcp_.time_prev[c]) * ge;
        }
        for (auto f : m.Faces()) {
          ffv_.iter_curr[f] +=
              (ffv_.time_curr[f] - 
               ffv_.time_prev[f]) * ge;
        }
      }
    }
  }
  // TODO: rewrite norm() using dist() where needed
  void MakeIteration() override {
    auto sem = m.GetSem("fluid-iter");
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
      ffv_.iter_prev = ffv_.iter_curr;

      CalcExtForce();

      CalcKinematicViscosity();

      timer_->Push("fluid.0.pressure-gradient");
      ffp_ = Interpolate(fcp_prev, mfcp_, m);
      fcgp_ = Gradient(ffp_, m);
      ffgp_ = Interpolate(fcgp_, mfcgp_, m);
      timer_->Pop();

      // initialize force for convdiff
      fc_cdf_.Reinit(m, Vect(0));
    }

    if (sem("explvisc")) {
      timer_->Push("fluid.1a.explicit-viscosity");
      // append viscous term
      for (size_t d = 0; d < dim; ++d) {
        fcwo_ = GetComponent(
            cd_->GetVelocity(Layers::iter_curr), d);
        auto ff = Interpolate(fcwo_, 
                              cd_->GetVelocityCond(d), m);
        auto gc = Gradient(ff, m);
        auto gf = Interpolate(gc, mfcf_, m); // adhoc: zero-der cond
        for (auto c : m.Cells()) {
          Vect s(0);
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            s += gf[f] * 
              (ffdk_[f] * m.GetOutwardSurface(c, q)[d]);
          }
          fc_cdf_[c] += s / m.GetVolume(c);
        }
      }
      timer_->Pop();
    }

    if (sem("forceappend")) {
      // append pressure gradient and external force
      for (auto c : m.AllCells()) {
        fc_cdf_[c] += fcgp_[c] * (-1.) + fcf_[c];
      }

      timer_->Push("fluid.2.convection-diffusion");
    }

    if (sem.Nested("convdiff-iter")) {
      cd_->MakeIteration();
    }

    if (sem("diag-comm")) {
      timer_->Pop();

      fck_.Reinit(m);
      for (auto c : m.Cells()) {
        Scal sum = 0.;
        for (size_t d = 0; d < dim; ++d) {
          sum += cd_->GetVelocityEquations(d)[c].CoeffSum();
        }
        fck_[c] = sum / dim;
      }

      m.Comm(&fck_);
    }
      
    if (sem("pcorr-assemble")) {
      // Define ffk_ on inner faces only
      ffk_ = Interpolate(fck_, mfck_, m);

      fcwe_ = cd_->GetVelocity(Layers::iter_curr);

      ffwe_ = Interpolate(
          fcwe_, mfcw_, m);

      // // needed for MMIM, now disabled
      //ff_velocity_iter_prev_ = Interpolate(
      //    cd_->GetVelocity(Layers::iter_prev),
      //    mfcw_, m);

      // Calc volumetric flux (asterisk)
      // using momentum interpolation (Rhie-Chow)
      // including hydrostatics-correction
      // TODO: Extend hydrostatics-correction on a non-uniform m
      timer_->Push("fluid.3.momentum-interpolation");
      ffve_.Reinit(m);
      const Scal rh = par->rhie;
      for (auto f : m.Faces()) {
        const auto volume_flux_interpolated =
            ffwe_[f].dot(m.GetSurface(f));
        if (!ffb_[f]) {
          IdxCell cm = m.GetNeighbourCell(f, 0);
          IdxCell cp = m.GetNeighbourCell(f, 1);
          Vect dm = m.GetVectToCell(f, 0);
          Vect dp = m.GetVectToCell(f, 1);
          const auto wide =
              (ffgp_[f] - fff_[f]).dot(m.GetSurface(f));
          const auto compact =
              ((fcp_prev[cp] - fcp_prev[cm]) /
              (dp - dm).norm() - (*fffp_)[f]) * m.GetArea(f);
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
      for (auto f : m.Faces()) {
        auto& expr = ffvc_[f];
        expr.Clear();
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        Vect dm = m.GetVectToCell(f, 0);
        Vect dp = m.GetVectToCell(f, 1);
        auto coeff = - m.GetArea(f) /
            ((dp - dm).norm() * ffk_[f]);
        if (ffb_[f]) {
          coeff = 0.;
        }
        expr.InsertTerm(-coeff, cm);
        expr.InsertTerm(coeff, cp);
        // adhoc for periodic
        expr.SortTerms(true);
        expr.SetConstant(ffve_[f]);
      }
      timer_->Pop();

      timer_->Push("fluid.5.pressure-system");
      for (auto c : m.Cells()) {
        auto& eqn = fcpcs_[c];
        Expr flux_sum;
        for (size_t i = 0; i < m.GetNumNeighbourFaces(c); ++i) {
          IdxFace f = m.GetNeighbourFace(c, i);
          flux_sum +=
              ffvc_[f] *
              m.GetOutwardFactor(c, i);
        }
        eqn =
            flux_sum -
            Expr((*fcsv_)[c] *
                 m.GetVolume(c));
      }

      // Account for cell conditions for pressure
      for (auto it = mccp_.cbegin();
          it != mccp_.cend(); ++it) {
        IdxCell c(it->GetIdx());
        ConditionCell* cond = it->GetValue().get();
        if (auto cond_value = dynamic_cast<ConditionCellValue<Scal>*>(cond)) {
          for (auto idxlocal : m.Cells()) {
            auto& eqn = fcpcs_[idxlocal];
            if (idxlocal == c) { 
              eqn.SetKnownValueDiag(c, cond_value->GetValue());
            } else {
              // Substitute value to obtain symmetrix matrix
              eqn.SetKnownValue(c, cond_value->GetValue());
            }
          }
        }
      }
    }

    if (sem("pcorr-solve")) {
      timer_->Pop();
      auto l = ConvertLs(fcpcs_, lsa_, lsb_, lsx_, m);
      using T = typename M::LS::T;
      l.t = T::symm;
      m.Solve(l);
      timer_->Push("fluid.6.pressure-solve");
    }

    if (sem("pcorr-comm")) {
      timer_->Pop();

      timer_->Push("fluid.7.correction");

      // Copy solution
      fcpc_.Reinit(m);
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
      for (auto c : m.Cells()) {
        fcp_curr[c] = fcp_prev[c] +
            pr * fcpc_[c];
      }
      m.Comm(&fcp_curr);

      fcgpc_ = Gradient(
          Interpolate(fcpc_, mfcpc_, m),
          m);

      // Correct the velocity
      fcvc_.Reinit(m);
      auto& u = cd_->GetVelocity(Layers::iter_curr);
      for (auto c : m.Cells()) {
        fcvc_[c] = fcgpc_[c] / (-fck_[c]);
      }
    }

    if (sem.Nested("convdiff-corr")) {
      // correct and comm
      cd_->CorrectVelocity(Layers::iter_curr, fcvc_);
    }

    if (sem("pcorr-fluxes")) {
      // Calc divergence-free volume fluxes
      for (auto f : m.Faces()) {
        ffv_.iter_curr[f] =
            ffvc_[f].Evaluate(fcpc_);
      }
      timer_->Pop();

      // TODO: SIMPLER removed

      this->IncIter();
      // m.Comm(&fc_velocity_curr);  // Comm done by CorrectVelocity
    }
  }
  void FinishStep() override {
    auto sem = m.GetSem("fluid-finish");
    if (sem("inctime")) {
      fcp_.time_prev = fcp_.time_curr;
      ffv_.time_prev = ffv_.time_curr;
      fcp_.time_curr = fcp_.iter_curr;
      ffv_.time_curr = ffv_.iter_curr;
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
    return ffv_.time_curr;
  }
  const FieldFace<Scal>& GetVolumeFlux(Layers layer) override {
    return ffv_.Get(layer);
  }
  double GetAutoTimeStep() override { 
    double dt = 1e10;
    auto& flux = ffv_.time_curr;
    for (auto c : m.Cells()) {
      for (size_t i = 0; i < m.GetNumNeighbourFaces(c); ++i) {
        IdxFace f = m.GetNeighbourFace(c, i);
        if (flux[f] != 0.) {
          dt = std::min<Scal>(
              dt, std::abs(m.GetVolume(c) / flux[f]));
        }
      }
    }
    return dt; 
  }
};

template <class M>
void FluidSimple<M>::Update(typename CD::Par& d, const Par& p) {
  // Update convdiff parameters
  d.relax = p.vrelax;
  d.guessextra = p.guessextra;
  d.second = p.second;
}

} // namespace solver

