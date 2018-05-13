#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "util/metrics.h"
#include "convdiffvi.h"
#include "fluid.h"

namespace solver {

// Rules:
// - Each function assumes that all fields on layers
//     iter_prev, time_prev, time_curr
//   and force, source and viscosity
//   are known and doesn't modify them. 
// - No function except for MakeIteration refers to iter_curr.

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
  using P::fcf_; // [i]
  using P::ffbp_; // [i]
  using P::fcsv_;
  using P::fcsm_;

  M& m; // mesh
  GRange<size_t> dr_;  // dimension range

  LayersData<FieldFace<Scal>> ffv_; // volume flux
  LayersData<FieldCell<Scal>> fcp_; // pressure

  std::shared_ptr<CD> cd_;

  // TODO: Const specifier for ConditionFace*
  
  // Face conditions
  MapFace<std::shared_ptr<ConditionFaceFluid>> mfc_; // fluid cond
  MapFace<std::shared_ptr<ConditionFace>> mfcw_; // velocity cond
  MapFace<std::shared_ptr<ConditionFace>> mfcp_; // pressure cond
  MapFace<std::shared_ptr<ConditionFace>> mfcf_; // force cond
  MapFace<std::shared_ptr<ConditionFace>> mfck_; // diag coeff cond
  MapFace<std::shared_ptr<ConditionFace>> mfcpc_; // pressure corr cond
  MapFace<std::shared_ptr<ConditionFace>> mfcd_; // dynamic viscosity cond

  // Cell conditions
  MapCell<std::shared_ptr<ConditionCellFluid>> mcc_; // fluid cell cond
  MapCell<std::shared_ptr<ConditionCell>> mccp_; // pressure cell cond
  MapCell<std::shared_ptr<ConditionCell>> mccw_; // velocity cell cond

  FieldFace<bool> ffbd_; // is boundary

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
  FieldCell<Scal> fcwo_;   // one velocity component
  FieldCell<Vect> fcb_;    // restored balanced force [s]
  FieldCell<Vect> fcfcd_;  // force for convdiff [i]

  // tmp
  std::array<FieldCell<Scal>, 3> fcta_; // TODO renname to vfct
  FieldCell<Scal> fct_;
  FieldCell<Scal> fct1_;
  FieldCell<Vect> fctv_;
  FieldFace<Vect> fftv_;
  FieldFace<Vect> fftv2_;
  FieldFace<Scal> fft_;

  // Face fields:
  FieldFace<Scal> ffd_;    // dynamic viscosity
  FieldFace<Scal> ffp_;    // pressure
  FieldFace<Vect> ffwe_;   // predicted velocity 
  FieldFace<Scal> ffve_;   // predicted volume flux [i]
  FieldFace<Scal> ffk_;    // diag coeff of velocity equation 
  FieldFace<Expr> ffvc_;   // expression for corrected volume flux [i]
  FieldFace<Vect> ffb_;    // restored balanced force [i]

  // LS
  std::vector<Scal> lsa_, lsb_, lsx_;

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

      if (auto cd = dynamic_cast<GivenPressure<M>*>(cb)) {
        auto pb = mccp_[i].get();
        auto pd = dynamic_cast<ConditionCellValueFixed<Scal>*>(pb);
        pd->Set(cd->GetPressure());
      } else if (auto cd = dynamic_cast<GivenVelocityAndPressure<M>*>(cb)) {
        auto pb = mccp_[i].get();
        auto pd = dynamic_cast<ConditionCellValueFixed<Scal>*>(pb);
        pd->Set(cd->GetPressure());
        auto wb = mccw_[i].get();
        auto wd = dynamic_cast<ConditionCellValueFixed<Vect>*>(wb);
        wd->Set(cd->GetVelocity());
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

  // Restore force from projections.
  // ffbp: force projections, bp=b.dot(n)
  // Output:
  // fcb, ffb: restored force
  // fcta_: modified tmp field
  void CalcExtForce(const FieldFace<Scal>& ffbp, 
                    FieldCell<Vect>& fcb, FieldFace<Vect>& ffb) {
    auto sem = m.GetSem("extforce");

    if (sem("loc")) {
      // XXX specific for Cartesian mesh
      // TODO consider just a weighted average of fn * n
      //      weight should be proportional to accuracy gradient approx
      //      which is better if surface area is larger
      fcb.Reinit(m);
      for (auto c : m.Cells()) {
        Vect s(0);
        for (auto q : m.Nci(c)) {
          // TODO: revise for non-rectangular cell
          IdxFace f = m.GetNeighbourFace(c, q);
          s += m.GetSurface(f) *
              (ffbp[f] * m.GetCenter(c).dist(m.GetCenter(f)));
        }
        fcb[c] = s / m.GetVolume(c);
      }

      // TODO: comm fffp_ instead
      for (auto d : dr_) {
        fcta_[d] = GetComponent(fcb, d);
        m.Comm(&fcta_[d]);
      }
    }

    if (sem("copy")) {
      for (auto d : dr_) {
        SetComponent(fcb, d, fcta_[d]);
      }
      ffb = Interpolate(fcb, mfcf_, m);
    }
  }

 public:
  struct Par {
    Scal vrelax = 0.8;   // velocity relaxation factor [0,1]
    Scal prelax = 1.;   // pressure relaxation factor [0,1]
    Scal rhie = 1.;     // Rhie-Chow factor [0,1] (0 disable, 1 full)
    bool second = true; // second order in time
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
              FieldCell<Scal>* fcr,  // density
              FieldCell<Scal>* fcd,  // dynamic viscosity
              FieldCell<Vect>* fcf,  // force 
              FieldFace<Scal>* ffbp, // balanced force projections 
              FieldCell<Scal>* fcsv, // volume source
              FieldCell<Scal>* fcsm, // mass source
              double time, double time_step,
              std::shared_ptr<Par> par
              )
      : FluidSolver<M>(time, time_step, fcr, fcd, fcf, ffbp, fcsv, fcsm)
      , m(m), dr_(0, dim), mfc_(mfc) , mcc_(mcc)
      , ffvc_(m), fcpcs_(m), par(par)
  {
    using namespace fluid_condition;

    ffbd_.Reinit(m, false);
    for (auto it : mfc_) {
      IdxFace i = it.GetIdx();
      ffbd_[i] = true;
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
              m, fcw, mfcw_, mccw_, fcr, &ffd_, 
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
      // rotate layers
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
  // Rhie-Chow interpolation of predicted volume flux
  // including balanced force (hydrostatics and surface tension)
  // fcw: predicted velocity field [s]
  // fcp: pressure field [s]
  // fcgp: gradient of pressure field [s]
  // fck, ffk: diag coeff [s]
  // Output:
  // ffv: result [i]
  // fftv_: modified tmp fields
  void RhieChow(const FieldCell<Vect>& fcw,
                const FieldCell<Scal>& fcp,  
                const FieldCell<Vect>& fcgp,
                const FieldCell<Scal>& fck,
                const FieldFace<Scal>& ffk,
                FieldFace<Scal>& ffv) {
    fftv_ = Interpolate(fcw, mfcw_, m); 

    const Scal rh = par->rhie; // rhie factor
    ffv.Reinit(m);
    for (auto f : m.Faces()) {
      // Init with mean flux
      ffv[f] = fftv_[f].dot(m.GetSurface(f));
      if (!ffbd_[f]) { // if not boundary
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        Vect dm = m.GetVectToCell(f, 0);
        Vect dp = m.GetVectToCell(f, 1);

        // compact pressure gradient
        Scal gp = (fcp[cp] - fcp[cm]) / (dp - dm).norm();

        // compact
        Scal o = ((*ffbp_)[f] - gp) * m.GetArea(f) / ffk[f];

        // wide
        Vect wm = (fcb_[cm] - fcgp_[cm]) / fck[cm];
        Vect wp = (fcb_[cp] - fcgp_[cp]) / fck[cp];
        Scal w = (wm + wp).dot(m.GetSurface(f)) * 0.5;

        // apply
        ffv[f] += rh * (o - w);
      } else { // if boundary
        // nop, keep mean flux
      }
    }

    // Apply meshvel
    for (auto f : m.Faces()) {
      ffv[f] -= par->meshvel.dot(m.GetSurface(f));
    }
  }
  // Apply cell conditions for pressure.
  // fcs: linear system in terms of correction of base pressure [i]
  // fcpb: base pressure [i]
  void ApplyPcCond(const FieldCell<Scal>& fcpb, FieldCell<Expr>& fcs) {
    for (auto it : mccp_) {
      IdxCell ct(it.GetIdx()); // cell target
      ConditionCell* cb = it.GetValue().get(); // cond base
      if (auto cd = dynamic_cast<ConditionCellValue<Scal>*>(cb)) {
        for (auto c : m.Cells()) {
          auto& e = fcs[c];
          Scal pc = cd->GetValue() - fcpb[c];
          if (c == ct) { 
            // Replace expression with [ct]-pc
            e.SetKnownValueDiag(ct, pc);
          } else {
            // Replace all ct terms with pc
            e.SetKnownValue(ct, pc);
          }
        }
      }
    }
  }
  // Flux expressions in terms of pressure correction pc:
  //   /  grad(pb+pc) * area / k + v, inner
  //   \  a, boundary
  // fcpb: base pressure [s]
  // ffk: diag coeff [i]
  // ffv: addition to flux [i]
  // Output:
  // ffe: result [i]
  void GetFlux(const FieldCell<Scal>& fcpb, const FieldFace<Scal>& ffk,
               const FieldFace<Scal>& ffv, FieldFace<Expr>& ffe) {
    ffe.Reinit(m);
    for (auto f : m.Faces()) {
      auto& e = ffe[f];
      e.Clear();
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      Vect dm = m.GetVectToCell(f, 0);
      Vect dp = m.GetVectToCell(f, 1);
      if (!ffbd_[f]) {  // inner
        Scal a = -m.GetArea(f) / ((dp - dm).norm() * ffk[f]);
        e.InsertTerm(-a, cm);
        e.InsertTerm(a, cp);
        e.SetConstant((fcpb[cp] - fcpb[cm]) * a + ffv[f]);
      } else { // boundary
        e.InsertTerm(0, cm);
        e.InsertTerm(0, cp);
        e.SetConstant(ffv[f]);
      }
    }
  }
  // Flux expressions in terms of pressure correction pc:
  //   /  grad(pc) * area / k + v, inner
  //   \  a, boundary
  // ffk: diag coeff [i]
  // ffv: addition to flux [i]
  // Output:
  // ffe: result [i]
  void GetFlux(const FieldFace<Scal>& ffk, const FieldFace<Scal>& ffv,
               FieldFace<Expr>& ffe) {
    ffe.Reinit(m);
    for (auto f : m.Faces()) {
      auto& e = ffe[f];
      e.Clear();
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      Vect dm = m.GetVectToCell(f, 0);
      Vect dp = m.GetVectToCell(f, 1);
      if (!ffbd_[f]) {  // inner
        Scal a = -m.GetArea(f) / ((dp - dm).norm() * ffk[f]);
        e.InsertTerm(-a, cm);
        e.InsertTerm(a, cp);
      } else { // boundary
        e.InsertTerm(0, cm);
        e.InsertTerm(0, cp);
      }
      e.SetConstant(ffv[f]);
    }
  }
  // Expressions for sum of fluxes and source:
  //   sum(v) - sv * vol
  // ffv: fluxes [i]
  // fcsv: volume source [i]
  // Output:
  // fce: result [i]
  void GetFluxSum(const FieldFace<Expr>& ffv, const FieldCell<Scal>& fcsv,
                  FieldCell<Expr>& fce) {
    fce.Reinit(m);
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      e.Clear();
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetNeighbourFace(c, q);
        e += ffv[f] * m.GetOutwardFactor(c, q);
      }
      e -= Expr(fcsv[c] * m.GetVolume(c));
    }
  }
  // Expressions for sum of fluxes.
  //   sum(v) 
  // ffv: fluxes [i]
  // Output:
  // fce: result [i]
  void GetFluxSum(const FieldFace<Expr>& ffv, FieldCell<Expr>& fce) {
    fce.Reinit(m);
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      e.Clear();
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetNeighbourFace(c, q);
        e += ffv[f] * m.GetOutwardFactor(c, q);
      }
    }
  }
  // Solve linear system fce = 0
  // fce: expressions [i]
  // Output:
  // fc: result [a]
  // lsa_, lsb_, lsx_: modified tmp fields
  void Solve(const FieldCell<Expr>& fce, FieldCell<Scal>& fc) {
    auto sem = m.GetSem("solve");
    if (sem("solve")) {
      auto l = ConvertLs(fcpcs_, lsa_, lsb_, lsx_, m);
      using T = typename M::LS::T; 
      l.t = T::symm; // solver type
      m.Solve(l);
    }
    if (sem("copy")) {
      fc.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc[c] = lsx_[i++];
      }
      m.Comm(&fc);
    }
  }
  // Get diagcoeff from current convdiff equations
  void GetDiagCoeff(FieldCell<Scal>& fck, FieldFace<Scal>& ffk) {
    auto sem = m.GetSem("diag");
    if (sem("local")) {
      fck.Reinit(m);
      for (auto c : m.Cells()) {
        Scal s = 0.;
        for (auto d : dr_) {
          // TODO consider separating diag from other coeffs
          //      use total sum for SIMPLEC only
          //sum += cd_->GetVelocityEquations(d)[c].CoeffSum();
          s += cd_->GetVelocityEquations(d)[c].Coeff(c); // XXX
        }
        fck[c] = s / dim;
      }

      m.Comm(&fck);
    }
    if (sem("interp")) {
      // TODO: remove mfck_ as probably not needed
      //       and add InterpolateInner
      ffk = Interpolate(fck, mfck_, m);
    }
  }
  // Append explicit part of viscous force.
  // fcw: velocity [a]
  // Output:
  // fcf += viscous term [i]
  // fftv_, fft_, fctv_, fftv2_: modified tmp fields
  void AppendExplViscous(const FieldCell<Vect>& fcw, FieldCell<Vect>& fcf) {
    auto& wc = fcw; // [a]
    auto& wf = fftv_; // [s]
    auto& wfo = fft_; // [s]
    auto& gc = fctv_; // [s]
    auto& gf = fftv2_; // [i]

    wf = Interpolate(wc, mfcw_, m); 
    for (auto d : dr_) {
      wfo = GetComponent(wf, d);
      gc = Gradient(wfo, m);
      gf = Interpolate(gc, mfcf_, m); // XXX adhoc zero-deriv cond
      for (auto c : m.Cells()) {
        Vect s(0);
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          s += gf[f] * (ffd_[f] * m.GetOutwardSurface(c, q)[d]);
        }
        fcf[c] += s / m.GetVolume(c);
      }
    }
  }
  // Restore pressure given velocity and volume flux
  // Assume MakeIteration() was called for convdiff solver
  // fcw: given velocity
  // ffv: given volume flux
  // fcpp: previous pressure
  // fcp: output pressure (may be aliased with fcpp)
  void CalcPressure(const FieldCell<Vect>& fcw,
                    const FieldFace<Scal>& ffv,
                    const FieldCell<Scal>& fcpp,
                    FieldCell<Scal>& fcp) {
    auto sem = m.GetSem("calcpressure");
    auto& fce = fctv_; // evaluation of velocity equations
    auto& ffa = fft_; // addition to flux TODO revise comment
    if (sem.Nested("cd-asm")) {
      cd_->Assemble(fcw, ffv);
    }
    if (sem.Nested("diag")) {
      GetDiagCoeff(fck_, ffk_);
    }
    if (sem("eval")) {
      for (auto d : dr_) {
        auto& w = fct_; // w^{s+1}
        w = GetComponent(fcw, d);
        auto& cd = cd_->GetSolver(d);
        auto& fce = cd.GetEquations();
        fcta_[d].Reinit(m);
        for (auto c : m.Cells()) {
          auto& e = fce[c];
          fcta_[d][c] = e.GetConstant(); // Evaluate(0)
        }
        m.Comm(&fcta_[d]);
      }
    }
    if (sem("assemble")) {
      // copy eval to vect
      fce.Reinit(m);
      for (auto d : dr_) {
        SetComponent(fce, d, fcta_[d]);
      }

      const Scal rh = par->rhie; // rhie factor
      ffa.Reinit(m);
      for (auto f : m.Faces()) {
        auto& e = ffvc_[f];
        e.Clear();

        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);

        if (!ffbd_[f]) { // if not boundary
          auto s = m.GetSurface(f);
          auto sa = m.GetArea(f);
          auto kf = rh * sa / ffk_[f];
          Vect bm = fcw[cm] - (fce[cm] - fcgp_[cm] + fcb_[cm]) / fck_[cm] * rh;
          Vect bp = fcw[cp] - (fce[cp] - fcgp_[cp] + fcb_[cp]) / fck_[cp] * rh;
          ffa[f] = (bm + bp).dot(s) * 0.5 + (*ffbp_)[f] * kf - ffv[f];
        } else { // if boundary
          ffa[f] = 0.;
        }
      }

      GetFlux(fcpp, ffk_, ffa, ffvc_);

      GetFluxSum(ffvc_, fcpcs_);

      ApplyPcCond(fcpp, fcpcs_);
    }

    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs_, fcpc_);
    }

    if (sem("apply")) {
      fcgpc_ = Gradient(Interpolate(fcpc_, mfcpc_, m), m);

      // Correct pressure
      for (auto c : m.Cells()) {
        fcp[c] = fcpp[c] + fcpc_[c]; 
      }
      m.Comm(&fcp);
    }
  }
  void UpdateBc(typename M::Sem& sem) {
    if (sem.Nested("bc-inletflux")) {
      UpdateInletFlux();
    }
    if (sem.Nested("bc-outlet")) {
      UpdateOutletBaseConditions();
    }
    if (sem("bc-derived")) {
      UpdateDerivedConditions();
    }
  }
  void CalcForce(typename M::Sem& sem) {
    if (sem.Nested("forceinit")) {
      CalcExtForce(*ffbp_, fcb_, ffb_);
    }

    if (sem("forceinit")) {
      // initialize force for convdiff
      fcfcd_.Reinit(m, Vect(0));
      // append force and balanced force
      for (auto c : m.Cells()) {
        fcfcd_[c] += (*fcf_)[c] + fcb_[c];
      }
    }

    if (sem("explvisc")) {
      AppendExplViscous(cd_->GetVelocity(Layers::iter_curr), fcfcd_);
    }
  }
  // TODO: rewrite norm() using dist() where needed
  void MakeIteration() override {
    auto sem = m.GetSem("fluid-iter");
    auto& fcp_prev = fcp_.iter_prev;
    auto& fcp_curr = fcp_.iter_curr;
    if (sem("init")) {
      // update convdiff par
      Update(*cd_->GetPar(), *par);
      // interpolate visosity 
      ffd_ = Interpolate(*fcd_, mfcd_, m, par->forcegeom);

      // rotate layers
      fcp_prev = fcp_curr;
      ffv_.iter_prev = ffv_.iter_curr;
    }

    UpdateBc(sem);

    CalcForce(sem);

    if (sem("pgrad")) {
      ffp_ = Interpolate(fcp_curr, mfcp_, m);
      fcgp_ = Gradient(ffp_, m);

      // append pressure gradient to force
      for (auto c : m.Cells()) {
        fcfcd_[c] += fcgp_[c] * (-1.);
      }
    }

    if (par->simpler) {
      if (sem.Nested("simpler")) {
        CalcPressure(cd_->GetVelocity(Layers::iter_curr), 
                     ffv_.iter_curr, fcp_curr, fcp_curr);
      }

      if (sem("pgrad")) {
        ffp_ = Interpolate(fcp_curr, mfcp_, m);
        fcgp_ = Gradient(ffp_, m);

        // append pressure correction gradient to force
        for (auto c : m.Cells()) {
          fcfcd_[c] += fcgpc_[c] * (-1.);
        }
      }
    }

    if (sem.Nested("convdiff-iter")) {
      // Solve for predictor velocity
      cd_->MakeIteration();
    }

    if (sem.Nested("diag")) {
      GetDiagCoeff(fck_, ffk_);
    }

    if (sem("pcorr-assemble")) {
      RhieChow(cd_->GetVelocity(Layers::iter_curr), 
               fcp_curr, fcgp_, fck_, ffk_, ffve_);

      GetFlux(ffk_, ffve_, ffvc_);

      GetFluxSum(ffvc_, *fcsv_, fcpcs_);

      ApplyPcCond(fcp_curr, fcpcs_);
    }

    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs_, fcpc_);
    }

    if (sem("pcorr-apply")) {
      if (!par->simpler) {
        // Correct pressure
        Scal pr = par->prelax; // pressure relaxation
        for (auto c : m.AllCells()) {
          fcp_curr[c] += pr * fcpc_[c];
        }
      }

      // Calc divergence-free volume fluxes
      for (auto f : m.Faces()) {
        ffv_.iter_curr[f] = ffvc_[f].Evaluate(fcpc_);
      }

      fcgpc_ = Gradient(Interpolate(fcpc_, mfcpc_, m), m);

      // Calc velocity correction
      fcwc_.Reinit(m);
      for (auto c : m.Cells()) {
        fcwc_[c] = fcgpc_[c] / (-fck_[c]);
      }
    }

    if (sem.Nested("convdiff-corr")) {
      // Correct velocity and comm
      cd_->CorrectVelocity(Layers::iter_curr, fcwc_); 
      //const_cast<FieldCell<Vect>&>(cd_->GetVelocity(Layers::iter_curr)) =
      //    cd_->GetVelocity(Layers::iter_prev); // XXX: adhoc keep velocity
    }

    if (sem("inc-iter")) {
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

