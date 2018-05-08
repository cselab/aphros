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

  // domain (cells/faces)
  // [i]: inner
  // [s]: support
  // [a]: all 
  

  using P::fcr_;
  using P::fcd_;
  using P::fffp_; // [i]
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

  // notation:
  // p: pressure
  // gp: pressure gradient
  // w: velocity
  // v: volume flux
  // we: predicted velocity (after solving velocity equations)
  // ve: predicted volume flux
  
  // Cell fields:
  FieldCell<Vect> fcgp_;   // gradient of pressure 
  FieldCell<Vect> fcwe_;   // predicted velocity 
  FieldCell<Scal> fck_;    // diag coeff of velocity equation 
  FieldCell<Expr> fcpcs_;  // pressure correction linear system [i]
  FieldCell<Scal> fcpc_;   // pressure correction
  FieldCell<Vect> fcgpc_;  // gradient of pressure correction
  FieldCell<Vect> fcwc_;   // velocity correction
  FieldCell<Scal> fcwo_;   // one velocitt component
  FieldCell<Scal> fcdk_;   // kinematic viscosity
  FieldCell<Vect> fcf_;    // restored vector force [s]
  FieldCell<Vect> fcfcd_;  // force for convdiff [i]

  FieldCell<Scal> fct0_;  // tmp
  FieldCell<Scal> fct1_;
  FieldCell<Scal> fct2_;

  // Face fields:
  FieldFace<Scal> ffp_;    // pressure
  FieldFace<Vect> ffgp_;   // gradient of pressure 
  FieldFace<Vect> ffwe_;   // predicted velocity 
  FieldFace<Scal> ffve_;   // predicted volume flux [i]
  FieldFace<Scal> ffk_;    // diag coeff of velocity equation 
  FieldFace<Expr> ffvc_;   // expression for corrected volume flux [i]
  FieldFace<Vect> fff_;    // restored vector force [i]
  FieldFace<Scal> ffdk_;   // kinematic viscosity

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
    auto sem = m.GetSem("extforce");

    if (sem("loc")) {
      timer_->Push("fluid.1.force-correction");
      // Restore vector force from faces
      fcf_.Reinit(m);
      for (auto c : m.Cells()) {
        Vect s(0);
        for (auto q : m.Nci(c)) {
          // TODO: revise for non-rectangular cell
          IdxFace f = m.GetNeighbourFace(c, q);
          s += m.GetSurface(f) *
              ((*fffp_)[f] * m.GetCenter(c).dist(m.GetCenter(f)));
        }
        fcf_[c] = s / m.GetVolume(c);
      }

      // TODO: comm fffp_ instead
      fct0_ = GetComponent(fcf_, 0);
      fct1_ = GetComponent(fcf_, 1);
      fct2_ = GetComponent(fcf_, 2);
      m.Comm(&fct0_);
      m.Comm(&fct1_);
      m.Comm(&fct2_);
    }

    if (sem("copy")) {
      SetComponent(fcf_, 0, fct0_);
      SetComponent(fcf_, 1, fct1_);
      SetComponent(fcf_, 2, fct2_);

      // Interpolate vector force to faces
      fff_ = Interpolate(fcf_, mfcf_, m);
      timer_->Pop();
    }
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
        mfcw_[i] = std::make_shared<
            ConditionFaceValueFixed<Vect>>(cd->GetVelocity(), nci);
        mfcp_[i] = std::make_shared<
            ConditionFaceExtrapolation>(nci);
      } else if (auto cd = dynamic_cast<Inlet<M>*>(cb)) {
        mfcw_[i] = std::make_shared<
            ConditionFaceValueFixed<Vect>>(cd->GetVelocity(), nci);
        mfcp_[i] = std::make_shared<
            ConditionFaceExtrapolation>(nci);
      } else if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
        mfcw_[i] = std::make_shared<
            ConditionFaceValueFixed<Vect>>(cd->GetVelocity(), nci);
        mfcp_[i] = std::make_shared<
            ConditionFaceExtrapolation>(nci);
      } else {
        throw std::runtime_error("Unknown fluid condition");
      }

      mfcgp_[i] = std::make_shared<
          ConditionFaceDerivativeFixed<Vect>>(Vect(0), nci);
      mfcf_[i] = std::make_shared<
          ConditionFaceDerivativeFixed<Vect>>(Vect(0), nci);
      mfcpc_[i] = std::make_shared<
          ConditionFaceExtrapolation>(nci);
      mfcd_[i] = std::make_shared<
          ConditionFaceDerivativeFixed<Scal>>(0., nci);
      mfck_[i] = std::make_shared<
          ConditionFaceDerivativeFixed<Scal>>(0, nci);
    }

    for (auto it : mcc_) {
      IdxCell c = it.GetIdx();
      ConditionCellFluid* cb = it.GetValue().get(); // cond base

      if (auto cd = dynamic_cast<GivenPressure<M>*>(cb)) {
        mccp_[c] = std::make_shared<
            ConditionCellValueFixed<Scal>>(cd->GetPressure());
      } else if (auto cd = dynamic_cast<GivenVelocityAndPressure<M>*>(cb)) {
        mccp_[c] = std::make_shared<
            ConditionCellValueFixed<Scal>>(cd->GetPressure());
        mccw_[c] = std::make_shared<
            ConditionCellValueFixed<Vect>>(cd->GetVelocity());
      } else {
        throw std::runtime_error("Unknown fluid cell condition");
      }
    }

    // Init convdiff solver
    {
      auto p = std::make_shared<typename CD::Par>();
      Update(*p, *par); // p from par

      fcfcd_.Reinit(m, Vect(0));
      cd_ = std::make_shared<
          ConvectionDiffusionImplicit<M>>(
              m, fcw, mfcw_, mccw_, fcr, &ffdk_, 
              &fcfcd_, &ffv_.iter_prev, time, time_step, p);
    }

    fcp_.time_curr.Reinit(m, 0.);
    fcp_.time_prev = fcp_.time_curr;

    // Calc initial volume fluxes
    fcwe_ = cd_->GetVelocity();
    ffwe_ = Interpolate(fcwe_, mfcw_, m);
    ffv_.time_curr.Reinit(m, 0.);
    for (auto f : m.Faces()) {
      ffv_.time_curr[f] = ffwe_[f].dot(m.GetSurface(f));
    }
    // Apply meshvel
    const Vect& meshvel = par->meshvel;
    for (auto f : m.Faces()) {
      ffv_.time_curr[f] -= meshvel.dot(m.GetSurface(f));
    }

    ffv_.time_prev = ffv_.time_curr;

    ffp_ = Interpolate(fcp_.time_curr, mfcp_, m);
    fcgp_ = Gradient(ffp_, m);
    ffgp_ = Interpolate(fcgp_, mfcgp_, m);
  }
  void StartStep() override {
    auto sem = m.GetSem("fluid-start");
    if (sem("convdiff-init")) {
      this->ClearIter();
      if (IsNan(fcp_.time_curr)) {
        throw std::runtime_error("simple::StartStep(): NaN pressure");
      }
      cd_->SetTimeStep(this->GetTimeStep());
    }

    if (sem.Nested("convdiff-start")) {
      cd_->StartStep();
    }

    if (sem("convdiff-start")) {
      fcp_.iter_curr = fcp_.time_curr;
      ffv_.iter_curr = ffv_.time_curr;
      // initial guess from extrapolation
      const Scal ge = par->guessextra;
      if (ge != 0.) {
        for (auto c : m.SuCells()) {
          fcp_.iter_curr[c] += (fcp_.time_curr[c] - fcp_.time_prev[c]) * ge;
        }
        for (auto f : m.Faces()) {
          ffv_.iter_curr[f] += (ffv_.time_curr[f] - ffv_.time_prev[f]) * ge;
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

    if (sem.Nested("extforce")) {
      CalcExtForce();
    }

    if (sem("pgrad")) {
      Update(*cd_->GetPar(), *par);
      UpdateDerivedConditions();

      fcp_prev = fcp_curr;
      ffv_.iter_prev = ffv_.iter_curr;

      CalcKinematicViscosity();

      timer_->Push("fluid.0.pressure-gradient");
      ffp_ = Interpolate(fcp_prev, mfcp_, m);
      fcgp_ = Gradient(ffp_, m);
      ffgp_ = Interpolate(fcgp_, mfcgp_, m);
      timer_->Pop();

      // initialize force for convdiff
      fcfcd_.Reinit(m, Vect(0));
    }

    if (sem("explvisc")) {
      timer_->Push("fluid.1a.explicit-viscosity");
      // append explicit part of viscous term
      for (size_t d = 0; d < dim; ++d) {
        fcwo_ = GetComponent(cd_->GetVelocity(Layers::iter_curr), d);
        auto ff = Interpolate(fcwo_, cd_->GetVelocityCond(d), m);
        auto gc = Gradient(ff, m);
        auto gf = Interpolate(gc, mfcf_, m); // adhoc: zero-der cond
        for (auto c : m.Cells()) {
          Vect s(0);
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            s += gf[f] * (ffdk_[f] * m.GetOutwardSurface(c, q)[d]);
          }
          fcfcd_[c] += s / m.GetVolume(c);
        }
      }
      timer_->Pop();
    }

    if (sem("forceappend")) {
      // append pressure gradient and external force
      for (auto c : m.Cells()) {
        fcfcd_[c] += fcgp_[c] * (-1.) + fcf_[c];
      }

      timer_->Push("fluid.2.convection-diffusion");
    }

    if (sem.Nested("convdiff-iter")) {
      // Solve for predictor velocity
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
      // TODO: remove mfck_ as probably not needed
      ffk_ = Interpolate(fck_, mfck_, m);

      fcwe_ = cd_->GetVelocity(Layers::iter_curr);

      ffwe_ = Interpolate(fcwe_, mfcw_, m);

      // Calc predicted volumetric flux 
      // using Momentum Interpolation (Rhie-Chow)
      // including balanced force (hydrostatics and surface tension)
      timer_->Push("fluid.3.momentum-interpolation");
      ffve_.Reinit(m);
      const Scal rh = par->rhie; // rhie factor
      for (auto f : m.Faces()) {
        // Init with interpolated flux
        ffve_[f] = ffwe_[f].dot(m.GetSurface(f));
        if (!ffb_[f]) { // if not boundary
          IdxCell cm = m.GetNeighbourCell(f, 0);
          IdxCell cp = m.GetNeighbourCell(f, 1);
          Vect dm = m.GetVectToCell(f, 0);
          Vect dp = m.GetVectToCell(f, 1);
          // TODO: rename correction to e.g. surplus 
    
          // Velocity equation structure:
          // w += f / k     // velocity correction
          // v += w * Surface  // volume flux correction
          // where f: force, k: diag coef (time step + density)

          // XXX: Consistency condition: 
          //      average of compact face gradients = cell gradient.

          // Wide approx of volume flux correction
          Scal qw = (fff_[f] - ffgp_[f]).dot(m.GetSurface(f));

          // Compact approx for pressure gradient
          Scal gp = (fcp_prev[cp] - fcp_prev[cm]) / (dp - dm).norm();
          // Compact approx of volume flux correction
          Scal qc = ((*fffp_)[f] - gp) * m.GetArea(f);

          // Corrected volume flux
          ffve_[f] += rh * (qc - qw) / ffk_[f];
        } else { // if boundary
          // nop, keep interpolated flux
        }
      }

      // Apply meshvel
      const Vect& meshvel = par->meshvel;
      for (auto f : m.Faces()) {
        ffve_[f] -= meshvel.dot(m.GetSurface(f));
      }

      timer_->Pop();

      timer_->Push("fluid.4.volume-flux");
      // Expressions for corrected volume flux
      // in terms of pressure
      // corrected = predicted + pressure gradient * area / diag coeff
      for (auto f : m.Faces()) {
        auto& e = ffvc_[f];
        e.Clear();
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        Vect dm = m.GetVectToCell(f, 0);
        Vect dp = m.GetVectToCell(f, 1);
        auto a = -m.GetArea(f) / ((dp - dm).norm() * ffk_[f]);
        if (ffb_[f]) { // keep on boundaries
          a = 0.;
        }
        e.InsertTerm(-a, cm);
        e.InsertTerm(a, cp);
        e.SetConstant(ffve_[f]);
      }
      timer_->Pop();

      timer_->Push("fluid.5.pressure-system");
      // System for pressure correction
      // sum of corrected volume fluxes + volume source = 0.
      for (auto c : m.Cells()) {
        auto& e = fcpcs_[c];
        Expr s;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          s += ffvc_[f] * m.GetOutwardFactor(c, q);
        }
        e = s - Expr((*fcsv_)[c] * m.GetVolume(c));
      }

      // Apply cell conditions for pressure
      // Traverse all expressions for every condition
      for (auto it : mccp_) {
        IdxCell cc(it.GetIdx()); // cell cond
        ConditionCell* cb = it.GetValue().get(); // cond base
        if (auto cd = dynamic_cast<ConditionCellValue<Scal>*>(cb)) {
          for (auto c : m.Cells()) {
            auto& e = fcpcs_[c];
            if (c == cc) { 
              // Replace expression with p[c] = cd->GetValue()
              e.SetKnownValueDiag(cc, cd->GetValue());
            } else {
              // Replace all cc terms with value
              e.SetKnownValue(cc, cd->GetValue());
            }
          }
        }
      }
    }

    if (sem("pcorr-solve")) {
      timer_->Pop();
      // Convert to LS format
      auto l = ConvertLs(fcpcs_, lsa_, lsb_, lsx_, m);
      using T = typename M::LS::T; 
      l.t = T::symm; // solver type
      // Solve system (add request)
      m.Solve(l);
      timer_->Push("fluid.6.pressure-solve");
    }

    if (sem("pcorr-comm")) {
      // System solved
      timer_->Pop();

      timer_->Push("fluid.7.correction");

      // Copy solution
      fcpc_.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fcpc_[c] = lsx_[i++];
      }
      
      // Comm pressure correction
      // (needed to compute gradients for flux correction)
      m.Comm(&fcpc_);
    }

    if (sem("pcorr-apply")) {
      // Correct pressure
      Scal pr = par->prelax; // pressure relaxation
      for (auto c : m.AllCells()) {
        fcp_curr[c] = fcp_prev[c] + pr * fcpc_[c];
      }
      m.Comm(&fcp_curr);

      fcgpc_ = Gradient(Interpolate(fcpc_, mfcpc_, m), m);

      // Compute velocity correction
      fcwc_.Reinit(m);
      for (auto c : m.Cells()) {
        fcwc_[c] = fcgpc_[c] / (-fck_[c]);
      }
    }

    if (sem.Nested("convdiff-corr")) {
      // Correct velocity and comm
      cd_->CorrectVelocity(Layers::iter_curr, fcwc_);
    }

    if (sem("pcorr-fluxes")) {
      // Calc divergence-free volume fluxes
      for (auto f : m.Faces()) {
        ffv_.iter_curr[f] = ffvc_[f].Evaluate(fcpc_);
      }
      timer_->Pop();

      // TODO: add SIMPLER or PISO

      this->IncIter();
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

