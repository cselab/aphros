#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "util/metrics.h"
#include "convdiffvi.h"
#include "convdiffi.h"
#include "convdiffe.h"
#include "fluid.h"
#include "debug/isnan.h"

#include <fstream>

// Rules:
// - Each function assumes that all fields on layers
//     iter_prev, time_prev, time_curr
//   and force, source and viscosity
//   are known and doesn't modify them. 
// - No function except for MakeIteration refers to iter_curr.
//
// domain (cells/faces)
// [i]: inner
// [s]: support
// [a]: all 
//
// notation:
// p: pressure
// gp: pressure gradient
// w: velocity
// v: volume flux
// we: predicted velocity (after solving velocity equations)
// ve: predicted volume flux

namespace solver {

template <class M_>
struct Simple<M_>::Imp {
  using Owner = Simple<M_>;
  using CD = ConvDiffVect<M>; // convdiff solver
  // convdiff solver implicit
  using CDI = ConvDiffVectImp<M, ConvDiffScalImp<M>>; 
  // convdiff solver explicit
  using CDE = ConvDiffVectImp<M, ConvDiffScalExp<M>>; 

  Imp(Owner* owner, const FieldCell<Vect>& fcw,
      const MapFace<std::shared_ptr<CondFaceFluid>>& mfc, 
      const MapCell<std::shared_ptr<CondCellFluid>>& mcc,
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m), dr_(0, m.GetEdim())
      , drr_(m.GetEdim(), dim)
      , mfc_(mfc), mcc_(mcc), fcpcs_(m), ffvc_(m)
  {
    using namespace fluid_condition;

    ffbd_.Reinit(m, false);

    mfcw_ = GetVelCond(m, mfc_);
    for (auto it : mfc_) {
      IdxFace i = it.GetIdx();
      ffbd_[i] = true;
      CondFaceFluid* cb = it.GetValue().get();
      size_t nci = cb->GetNci();

      mfcf_[i]  = std::make_shared<CondFaceGradFixed<Vect>>(Vect(0), nci);
      mfcp_[i]  = std::make_shared<CondFaceExtrap>(nci);
      mfcpc_[i] = std::make_shared<CondFaceExtrap>(nci);
      mfcd_[i]  = std::make_shared<CondFaceGradFixed<Scal>>(0., nci);

      if (auto cd = dynamic_cast<NoSlipWall<M>*>(cb)) {
        // nop
      } else if (auto cd = dynamic_cast<Inlet<M>*>(cb)) {
        // nop
      } else if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
        // nop
      } else if (auto cd = dynamic_cast<SlipWall<M>*>(cb)) {
        mfcf_[i]  = std::make_shared<CondFaceReflect>(nci);
        mfcp_[i]  = std::make_shared<CondFaceGradFixed<Scal>>(0., nci);
        mfcpc_[i] = std::make_shared<CondFaceGradFixed<Scal>>(0, nci);
      } else {
        throw std::runtime_error("simple: unknown condition");
      }
    }

    for (auto it : mcc_) {
      IdxCell c = it.GetIdx();
      CondCellFluid* cb = it.GetValue().get(); // cond base

      if (auto cd = dynamic_cast<GivenPressure<M>*>(cb)) {
        mccp_[c] = std::make_shared<
            CondCellValFixed<Scal>>(cd->GetPressure());
      } else if (auto cd = dynamic_cast<GivenVelocityAndPressure<M>*>(cb)) {
        mccp_[c] = std::make_shared<
            CondCellValFixed<Scal>>(cd->GetPressure());
        mccw_[c] = std::make_shared<
            CondCellValFixed<Vect>>(cd->GetVelocity());
      } else {
        throw std::runtime_error("simple: unknown cell condition");
      }
    }

    InitConvDiff(fcw);

    fcp_.time_curr.Reinit(m, 0.);
    fcp_.time_prev = fcp_.time_curr;

    // Calc initial volume fluxes
    auto ffwe = Interpolate(cd_->GetVelocity(), mfcw_, m);
    ffv_.time_curr.Reinit(m, 0.);
    for (auto f : m.Faces()) {
      ffv_.time_curr[f] = ffwe[f].dot(m.GetSurface(f));
    }
    // Apply meshvel
    const Vect& meshvel = par->meshvel;
    for (auto f : m.Faces()) {
      ffv_.time_curr[f] -= meshvel.dot(m.GetSurface(f));
    }

    ffv_.time_prev = ffv_.time_curr;

    auto ffp = Interpolate(fcp_.time_curr, mfcp_, m);
    fcgp_ = Gradient(ffp, m);
  }

  void InitConvDiff(const FieldCell<Vect>& fcw) {
    fcfcd_.Reinit(m, Vect(0));

    if (par->conv == Conv::imp) {
      auto p = std::make_shared<typename CDI::Par>();
      Update(*p, *par);

      cd_ = std::make_shared<CDI>(
          m, fcw, mfcw_, mccw_, owner_->fcr_, &ffd_, &fcfcd_, 
          &ffv_.iter_prev, owner_->GetTime(), owner_->GetTimeStep(), p);
    } else if (par->conv == Conv::exp) {
      auto p = std::make_shared<typename CDE::Par>();
      Update(*p, *par);

      cd_ = std::make_shared<CDE>(
          m, fcw, mfcw_, mccw_, owner_->fcr_, &ffd_, &fcfcd_, 
          &ffv_.iter_prev, owner_->GetTime(), owner_->GetTimeStep(), p);
    }
  }

  // TODO: somehow track dependencies to define execution order
  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    // Face conditions
    for (auto it : mfc_) {
      IdxFace i = it.GetIdx();
      CondFaceFluid* cb = it.GetValue().get();
      auto p = mfcw_[i].get();

      if (auto cd = dynamic_cast<NoSlipWall<M>*>(cb)) {
        auto pd = dynamic_cast<CondFaceValFixed<Vect>*>(p);
        pd->Set(cd->GetVelocity());
      } else if (auto cd = dynamic_cast<Inlet<M>*>(cb)) {
        auto pd = dynamic_cast<CondFaceValFixed<Vect>*>(p);
        pd->Set(cd->GetVelocity());
      } else if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
        auto pd = dynamic_cast<CondFaceValFixed<Vect>*>(p);
        pd->Set(cd->GetVelocity());
      } else if (auto cd = dynamic_cast<SlipWall<M>*>(cb)) {
        // nop
      } else {
        throw std::runtime_error("simple: Unknown fluid condition");
      }
    }

    // Cell conditions
    for (auto it : mcc_) {
      IdxCell i = it.GetIdx();
      CondCellFluid* cb = it.GetValue().get();

      if (auto cd = dynamic_cast<GivenPressure<M>*>(cb)) {
        auto pb = mccp_[i].get();
        auto pd = dynamic_cast<CondCellValFixed<Scal>*>(pb);
        pd->Set(cd->GetPressure());
      } else if (auto cd = dynamic_cast<GivenVelocityAndPressure<M>*>(cb)) {
        auto pb = mccp_[i].get();
        auto pd = dynamic_cast<CondCellValFixed<Scal>*>(pb);
        pd->Set(cd->GetPressure());
        auto wb = mccw_[i].get();
        auto wd = dynamic_cast<CondCellValFixed<Vect>*>(wb);
        wd->Set(cd->GetVelocity());
      } else {
        throw std::runtime_error("simple: Unknown fluid cell condition");
      }
    }
  }
  /*
  // XXX adhoc, TODO: revise
  void UpdateWall() {
    using namespace fluid_condition;

    auto& vel = GetVelocity(Layers::iter_curr);
    auto sem = m.GetSem("");
    if (sem("")) {
      for (size_t id = 0; id < nid; ++id) {
        for (auto it : mfc_) {
          IdxFace i = it.GetIdx();
          CondFaceFluid* cb = it.GetValue().get(); // cond base

          if (auto cd = dynamic_cast<NoSlipWall<M>*>(cb)) {
            size_t nci = cd->GetNci();
            Scal w = (nci == 0 ? -1. : 1.);
            Vect n = m.GetNormal(i);
            cd->SetVelocity(cd->GetVelocity() + n * (dv * w));
          }
        }
      }
    }
  }
  */
  // TODO: consider updating from predictor velocity
  void UpdateInletFlux() {
    using namespace fluid_condition;
    size_t& nid = par->inletflux_numid;

    auto sem = m.GetSem("inletflux");
    if (sem("local")) { 
      ilft_.resize(nid);
      ilfe_.resize(nid);
      ila_.resize(nid);

      for (size_t id = 0; id < nid; ++id) {
        ilft_[id] = 0.;
        ilfe_[id] = 0.;
        ila_[id] = 0.;
      }

      // Extrapolate velocity to inlet from neighbour cells
      // and compute total fluxes
      auto& vel = GetVelocity(Layers::iter_curr);
      for (auto it : mfc_) {
        IdxFace i = it.GetIdx();
        CondFaceFluid* cb = it.GetValue().get(); // cond base

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
      
      for (size_t id = 0; id < nid; ++id) {
        m.Reduce(&ilft_[id], "sum");
        m.Reduce(&ilfe_[id], "sum");
        m.Reduce(&ila_[id], "sum");
      }
    }

    if (sem("corr")) {
      for (size_t id = 0; id < nid; ++id) {
        // Apply additive correction
        Scal dv = (ilft_[id] - ilfe_[id]) / ila_[id];  // velocity
        for (auto it : mfc_) {
          IdxFace i = it.GetIdx();
          CondFaceFluid* cb = it.GetValue().get(); // cond base

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
        CondFaceFluid* cb = it->GetValue().get(); // cond base

        size_t id = cb->GetNci();
        IdxCell c = m.GetNeighbourCell(i, id);
        if (m.IsInner(c)) {
          if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
            Scal w = (id == 0 ? 1. : -1.);
            Vect vc = GetVelocity(Layers::iter_curr)[c];
            Vect s = m.GetSurface(i);
            // clip normal component, let only positive
            vc -= s * (w * std::min(0., vc.dot(s) * w)  / s.dot(s));
            cd->SetVelocity(vc);
            fo += cd->GetVelocity().dot(s) * w;
            ao += m.GetArea(i);
          } else if (auto cd = dynamic_cast<Inlet<M>*>(cb)) {
            Scal w = (id == 0 ? -1. : 1.);
            fi += cd->GetVelocity().dot(m.GetSurface(i)) * w;
          }
        }
      }
      
      // Append volume source to inlet flux
      auto& fcsv = *owner_->fcsv_;
      for (auto i : m.Cells()) {
        fi += fcsv[i] * m.GetVolume(i);
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
        CondFaceFluid* cb = it->GetValue().get(); // cond base

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
  // fcb: restored force
  void CalcExtForce(const FieldFace<Scal>& ffbp, 
                    FieldCell<Vect>& fcb) {
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
              (ffbp[f] * m.GetVolume(c) / m.GetArea(f) * 0.5);
        }
        fcb[c] = s / m.GetVolume(c);
      }

      // TODO: comm ffbp_ instead
      m.Comm(&fcb);
    }
  }
  void Update(typename CDI::Par& d, const Par& p) {
    // Update convdiff parameters
    d.relax = p.vrelax;
    d.guessextra = p.guessextra;
    d.second = p.second;
    d.sc = p.convsc;
    d.df = p.convdf;
    d.linreport = p.linreport;
  }
  void Update(CDI* cd, const Par& p) {
    Update(*cd->GetPar(), p);
  }
  void Update(typename CDE::Par& d, const Par& p) {
    // Update convdiff parameters
    d.relax = p.vrelax;
    d.guessextra = p.guessextra;
    d.second = p.second;
    d.sc = p.convsc;
    d.df = p.convdf;
    d.linreport = p.linreport;
  }
  void Update(CDE* cd, const Par& p) {
    Update(*cd->GetPar(), p);
  }
  void StartStep() {
    auto sem = m.GetSem("fluid-start");
    if (sem("convdiff-init")) {
      owner_->ClearIter();
      CHECKNAN(fcp_.time_curr, m.CN())
      cd_->SetTimeStep(owner_->GetTimeStep());
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
  void RhieChow(const FieldCell<Vect>& fcw,
                const FieldCell<Scal>& fcp,  
                const FieldCell<Vect>& fcgp,
                const FieldCell<Scal>& fck,
                const FieldFace<Scal>& ffk,
                FieldFace<Scal>& ffv) {
    auto fftv = Interpolate(fcw, mfcw_, m); 

    const Scal rh = par->rhie; // rhie factor
    ffv.Reinit(m);
    auto& ffbp = *owner_->ffbp_;
    for (auto f : m.Faces()) {
      // Init with mean flux
      ffv[f] = fftv[f].dot(m.GetSurface(f));
      if (!ffbd_[f]) { // if not boundary
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);

        // compact pressure gradient
        Scal hr = m.GetArea(f) / m.GetVolume(cp);
        Scal gp = (fcp[cp] - fcp[cm]) * hr;

        // compact
        Scal o = (ffbp[f] - gp) * m.GetArea(f) / ffk[f];

        // wide
        Vect wm = (fcb_[cm] - fcgp[cm]) / fck[cm];
        Vect wp = (fcb_[cp] - fcgp[cp]) / fck[cp];
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
      CondCell* cb = it.GetValue().get(); // cond base
      // TODO: revise without loop over cells
      // TODO: adhoc, replced assuming no references from other cells
      //       (i.e. only excluded cells)
      if (auto cd = dynamic_cast<CondCellVal<Scal>*>(cb)) {
        auto& e = fcs[ct];
        Scal pc = cd->GetValue() - fcpb[ct];
        e.SetKnownValueDiag(ct, pc);
        /*
        for (auto c : m.Cells()) {
          auto& e = fcs[c];
          if (c == ct) { 
            // Replace expression with [ct]-pc
            e.SetKnownValueDiag(ct, pc);
          } else {
            // Replace all ct terms with pc
            e.SetKnownValue(ct, pc);
          }
        }
        */
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
      Scal hr = m.GetArea(f) / m.GetVolume(cp);
      if (!ffbd_[f]) {  // inner
        Scal a = -m.GetArea(f) * hr / ffk[f];
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
      Scal hr = m.GetArea(f) / m.GetVolume(cp);
      if (!ffbd_[f]) {  // inner
        Scal a = -m.GetArea(f) * hr / ffk[f];
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
  // m.GetSolveTmp(): modified temporary fields
  void Solve(const FieldCell<Expr>& fce, FieldCell<Scal>& fc) {
    auto sem = m.GetSem("solve");
    if (sem("solve")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);
      auto l = ConvertLs(fce, *lsa, *lsb, *lsx, m);

      {
        static int fi = 0;
        ++fi;
        std::ofstream f("lin_good_" + std::to_string(fi) + ".dat");

        using T = GVect<Scal, M::dim * 2 + 2>;
        FieldCell<T> fc(m, T(0));
        for (auto c : m.Cells()) {
          auto e = fce[c];
          size_t i = 0;
          T v;
          v[i++] = e[3].a;
          v[i++] = e[0].a;
          v[i++] = e[1].a;
          v[i++] = e[2].a;
          v[i++] = e[4].a;
          v[i++] = e[5].a;
          v[i++] = e[6].a;
          v[i++] = e.GetConstant();
          fc[c] = v;
        }

        f << fc;

        size_t j = 0;
        for (size_t i = 0; i < lsb->size(); ++i) {
          f << i << " "
            << (*lsa)[j + 3] << " "
            << (*lsa)[j + 0] << " "
            << (*lsa)[j + 1] << " "
            << (*lsa)[j + 2] << " "
            << (*lsa)[j + 4] << " "
            << (*lsa)[j + 5] << " "
            << (*lsa)[j + 6] << "\n "
            << (*lsb)[i] << " "
            << (*lsx)[i] << " "
            << "\n";
          j += 7;
        }

      }

      using T = typename M::LS::T; 
      l.t = T::symm; // solver type
      m.Solve(l);
    }
    if (sem("copy")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);

      {
        static int fi = 0;
        ++fi;
        std::ofstream f("lin_good_x_" + std::to_string(fi) + ".dat");

        for (size_t i = 0; i < lsx->size(); ++i) {
          f << i << " "
            << (*lsx)[i] << " "
            << "\n";
        }
      }

      fc.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc[c] = (*lsx)[i++];
      }
      CHECKNAN(fc, m.CN());
      m.Comm(&fc);
      if (par->linreport && m.IsRoot()) {
        std::cout 
            << "pcorr:"
            << " res=" << m.GetResidual()
            << " iter=" << m.GetIter()
            << std::endl;
      }
    }
  }
  // Get diagcoeff from current convdiff equations
  void GetDiagCoeff(FieldCell<Scal>& fck, FieldFace<Scal>& ffk) {
    auto sem = m.GetSem("diag");
    if (sem("local")) {
      fck.Reinit(m, 0);
      for (auto d : dr_) {
        auto fct = cd_->GetDiag(d);
        for (auto c : m.Cells()) {
          fck[c] += fct[c];
        }
      }
      for (auto c : m.Cells()) {
        fck[c] /= dr_.size();
      }

      CHECKNAN(fck, m.CN())

      m.Comm(&fck);
    }
    if (sem("interp")) {
      ffk.Reinit(m);
      InterpolateI(fck, ffk, m);
    }
  }
  // Append explicit part of viscous force.
  // fcw: velocity [a]
  // Output:
  // fcf += viscous term [i]
  void AppendExplViscous(const FieldCell<Vect>& fcw, FieldCell<Vect>& fcf) {
    auto wf = Interpolate(fcw, mfcw_, m); 
    for (auto d : dr_) {
      auto wfo = GetComponent(wf, d);
      auto gc = Gradient(wfo, m);
      auto gf = Interpolate(gc, mfcf_, m); // XXX adhoc zero-deriv cond
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
  // Output:
  // fcp: output pressure (may be aliased with fcpp)
  // fctv_: modified tmp 
  void CalcPressure(const FieldCell<Vect>& fcw,
                    const FieldFace<Scal>& ffv,
                    const FieldCell<Scal>& fcpp,
                    FieldCell<Scal>& fcp) {
    auto sem = m.GetSem("calcpressure");
    auto& fcl = fctv_; // evaluation of velocity equations
    if (sem.Nested("cd-asm")) {
      cd_->Assemble(fcw, ffv);
    }
    if (sem.Nested("diag")) {
      GetDiagCoeff(fck_, ffk_);
    }
    if (sem("eval")) {
      fcl.Reinit(m);
      for (auto d : dr_) {
        SetComponent(fcl, d, cd_->GetConst(d));
      }
      for (auto d : drr_) {
        SetComponent(fcl, d, 0);
      }
      m.Comm(&fcl);
    }
    if (sem("assemble")) {
      const Scal rh = par->rhie; // rhie factor
      FieldFace<Scal> ffa(m); // addition to flux TODO revise comment
      auto& ffbp = *owner_->ffbp_;
      for (auto f : m.Faces()) {
        auto& e = ffvc_[f];
        e.Clear();

        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);

        if (!ffbd_[f]) { // if not boundary
          auto s = m.GetSurface(f);
          auto sa = m.GetArea(f);
          auto kf = rh * sa / ffk_[f];
          Vect bm = fcw[cm] - (fcl[cm] - fcgp_[cm] + fcb_[cm]) / fck_[cm] * rh;
          Vect bp = fcw[cp] - (fcl[cp] - fcgp_[cp] + fcb_[cp]) / fck_[cp] * rh;
          ffa[f] = (bm + bp).dot(s) * 0.5 + ffbp[f] * kf - ffv[f];
        } else { // if boundary
          ffa[f] = 0.;
        }
      }
      fcl.Free();

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
      Scal pr = par->prelax; // pressure relaxation
      for (auto c : m.Cells()) {
        fcp[c] = fcpp[c] + fcpc_[c] * pr; 
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
      CalcExtForce(*owner_->ffbp_, fcb_);
    }

    if (sem("forceinit")) {
      // initialize force for convdiff
      fcfcd_.Reinit(m, Vect(0));
      // append force and balanced force
      auto& fcf = *owner_->fcf_;
      for (auto c : m.Cells()) {
        fcfcd_[c] += fcf[c] + fcb_[c];
      }
    }

    if (sem("explvisc")) {
      AppendExplViscous(cd_->GetVelocity(Layers::iter_curr), fcfcd_);
    }
  }
  // TODO: rewrite norm() using dist() where needed
  void MakeIteration() {
    auto sem = m.GetSem("fluid-iter");
    auto& fcp_prev = fcp_.iter_prev;
    auto& fcp_curr = fcp_.iter_curr;
    if (sem("init")) {
      // update convdiff par
      if (auto p = dynamic_cast<CDI*>(cd_.get())) {
        Update(p, *par);
      } else if (auto p = dynamic_cast<CDE*>(cd_.get())) {
        Update(p, *par);
      } else {
        throw std::runtime_error(sem.GetName() + ": unknown cd_");
      }
      // interpolate visosity 
      ffd_ = Interpolate(*owner_->fcd_, mfcd_, m);

      // rotate layers
      fcp_prev = fcp_curr;
      ffv_.iter_prev = ffv_.iter_curr;
    }

    UpdateBc(sem);

    CalcForce(sem);

    if (sem("pgrad")) {
      auto ffp = Interpolate(fcp_curr, mfcp_, m);
      fcgp_ = Gradient(ffp, m);

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
        auto ffp = Interpolate(fcp_curr, mfcp_, m);
        fcgp_ = Gradient(ffp, m);

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
      CHECKNAN(ffve_, m.CN())

      GetFlux(ffk_, ffve_, ffvc_);
      CHECKNAN(ffvc_, m.CN())

      GetFluxSum(ffvc_, *owner_->fcsv_, fcpcs_);
      CHECKNAN(fcpcs_, m.CN())

      ApplyPcCond(fcp_curr, fcpcs_);
    }

    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs_, fcpc_);
    }

    if (sem("pcorr-apply")) {
      CHECKNAN(fcpc_, m.CN())

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
      CHECKNAN(ffv_.iter_curr, m.CN())

      auto ffpc = Interpolate(fcpc_, mfcpc_, m);
      CHECKNAN(ffpc, m.CN())
      fcgpc_ = Gradient(ffpc, m);
      CHECKNAN(fcgpc_, m.CN())

      // Calc velocity correction
      fcwc_.Reinit(m);
      for (auto c : m.Cells()) {
        fcwc_[c] = fcgpc_[c] / (-fck_[c]);
      }
      CHECKNAN(fcwc_, m.CN())
    }

    if (sem.Nested("convdiff-corr")) {
      // Correct velocity and comm
      cd_->CorrectVelocity(Layers::iter_curr, fcwc_); 
    }

    if (sem("inc-iter")) {
      owner_->IncIter();
      fcwc_.Free();
      fck_.Free();
    }
  }
  void FinishStep() {
    auto sem = m.GetSem("fluid-finish");
    if (sem("inctime")) {
      fcp_.time_prev = fcp_.time_curr;
      ffv_.time_prev = ffv_.time_curr;
      fcp_.time_curr = fcp_.iter_curr;
      ffv_.time_curr = ffv_.iter_curr;
      CHECKNAN(fcp_.time_curr, m.CN())
      owner_->IncTime();
    }
    if (sem.Nested("convdiff-finish")) {
      cd_->FinishStep();
    }
  }
  double GetAutoTimeStep() { 
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
  const FieldCell<Vect>& GetVelocity(Layers l) const {
    return cd_->GetVelocity(l);
  }

  Owner* owner_;
  std::shared_ptr<Par> par;
  M& m; // mesh
  GRange<size_t> dr_;  // effective dimension range
  GRange<size_t> drr_;  // remaining dimensions

  // Face conditions
  MapFace<std::shared_ptr<CondFaceFluid>> mfc_; // fluid cond
  MapFace<std::shared_ptr<CondFace>> mfcw_; // velocity cond
  MapFace<std::shared_ptr<CondFace>> mfcp_; // pressure cond
  MapFace<std::shared_ptr<CondFace>> mfcf_; // force cond
  MapFace<std::shared_ptr<CondFace>> mfcpc_; // pressure corr cond
  MapFace<std::shared_ptr<CondFace>> mfcd_; // dynamic viscosity cond

  // Cell conditions
  MapCell<std::shared_ptr<CondCellFluid>> mcc_; // fluid cell cond
  MapCell<std::shared_ptr<CondCell>> mccp_; // pressure cell cond
  MapCell<std::shared_ptr<CondCell>> mccw_; // velocity cell cond

  LayersData<FieldFace<Scal>> ffv_; // volume flux
  LayersData<FieldCell<Scal>> fcp_; // pressure

  std::shared_ptr<CD> cd_;

  // TODO: Const specifier for CondFace*
  
  FieldFace<bool> ffbd_; // is boundary

  // used by UpdateOutletBaseConditions():
  Scal olfi_; // inlet flux
  Scal olfo_; // outlet flux
  Scal olao_; // outlet area
  // used by UpdateInletFlux():
  std::vector<Scal> ilft_; // target flux
  std::vector<Scal> ilfe_; // extrapolated flux
  std::vector<Scal> ila_; // area

  // Cell fields:
  FieldCell<Vect> fcgp_;   // gradient of pressure 
  FieldCell<Scal> fck_;    // diag coeff of velocity equation 
  FieldCell<Expr> fcpcs_;  // pressure correction linear system [i]
  FieldCell<Scal> fcpc_;   // pressure correction
  FieldCell<Vect> fcgpc_;  // gradient of pressure correction
  FieldCell<Vect> fcwc_;   // velocity correction
  FieldCell<Vect> fcb_;    // restored balanced force [s]
  FieldCell<Vect> fcfcd_;  // force for convdiff [i]

  // tmp
  FieldCell<Scal> fct_;
  FieldCell<Vect> fctv_;

  // Face fields:
  FieldFace<Scal> ffd_;    // dynamic viscosity
  FieldFace<Scal> ffve_;   // predicted volume flux [i]
  FieldFace<Scal> ffk_;    // diag coeff of velocity equation 
  FieldFace<Expr> ffvc_;   // expression for corrected volume flux [i]

};

template <class M_>
Simple<M_>::Simple(
    M& m, const FieldCell<Vect>& fcw,
    const MapFace<std::shared_ptr<CondFaceFluid>>& mfc,
    const MapCell<std::shared_ptr<CondCellFluid>>& mcc,
    FieldCell<Scal>* fcr, FieldCell<Scal>* fcd, 
    FieldCell<Vect>* fcf, FieldFace<Scal>* ffbp,
    FieldCell<Scal>* fcsv, FieldCell<Scal>* fcsm,
    double t, double dt, std::shared_ptr<Par> par)
    : FluidSolver<M>(t, dt, m, fcr, fcd, fcf, ffbp, fcsv, fcsm)
    , imp(new Imp(this, fcw, mfc, mcc, par))
{}

template <class M_>
Simple<M_>::~Simple() = default;

template <class M_>
auto Simple<M_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_>
void Simple<M_>::StartStep() {
  return imp->StartStep();
}

template <class M_>
void Simple<M_>::MakeIteration() {
  return imp->MakeIteration();
}

template <class M_>
void Simple<M_>::FinishStep() {
  return imp->FinishStep();
}

template <class M_>
auto Simple<M_>::GetVelocity(Layers l) const -> const FieldCell<Vect>& {
  return imp->GetVelocity(l);
}

template <class M_>
auto Simple<M_>::GetPressure(Layers l) const -> const FieldCell<Scal>& {
  return imp->fcp_.Get(l);
}

template <class M_>
auto Simple<M_>::GetVolumeFlux(Layers l) const -> const FieldFace<Scal>& {
  return imp->ffv_.Get(l);
}

template <class M_>
double Simple<M_>::GetAutoTimeStep() const {
  return imp->GetAutoTimeStep();
}

template <class M_>
double Simple<M_>::GetError() const {
  return imp->cd_->GetError();
}

template <class M_>
auto Simple<M_>::GetVelocityCond() const -> 
    const MapFace<std::shared_ptr<CondFace>>& {
  return imp->mfcw_;
}


} // namespace solver

