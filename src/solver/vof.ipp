#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>

#include "vof.h"
#include "geom/block.h"
#include "dump/vtk.h"
#include "reconst.h"
#include "partstr.h"
#include "normal.h"

namespace solver {

template <class M_>
struct Vof<M_>::Imp {
  using R = Reconst<Scal>;
  using PS = PartStr<Scal>;
  static constexpr size_t dim = M::dim;

  Imp(Vof* owner, const FieldCell<Scal>& fcu,
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), fca_(m, 0), fcn_(m, Vect(0)), fck_(m, 0), fckp_(m, 0)
  {
    fcu_.time_curr = fcu;
    for (auto it : mfc_) {
      IdxFace f = it.GetIdx();
      mfvz_[f] = std::make_shared<
          CondFaceGradFixed<Vect>>(Vect(0), it.GetValue()->GetNci());
    }

    // particle strings
    auto p = std::make_shared<typename PS::Par>();
    Update(p.get());
    partstr_ = std::unique_ptr<PS>(new PS(p));
  }
  void Update(typename PS::Par* p) const {
    Scal hc = GetCellSize().norminf(); // cell size

    p->leq = par->part_h;
    p->relax = par->part_relax;
    p->constr = par->part_constr;
    p->npmax = par->part_np;
    p->segcirc = par->part_segcirc;
    p->hc = hc;

    p->kstr = par->part_kstr;
    p->kbend = par->part_kbend;
    p->kattr = par->part_kattr;
    p->bendmean = par->part_bendmean;
    p->ka = par->part_kattr;
    p->kt = par->part_kbend;
    p->kx = par->part_kstr;
    p->tmax = par->part_tmax;
    p->dtmax = par->part_dtmax;
    p->anglim = par->part_anglim;
  }
  // Dump cut polygons
  void DumpPoly() {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      auto h = GetCellSize();
      for (auto c : m.Cells()) {
        const Scal th = par->poly_intth;
        Scal u = fcu_.iter_curr[c];
        if (fci_[c] && u > th && u < 1. - th) {
          dl_.push_back(R::GetCutPoly(m.GetCenter(c), fcn_[c], fca_[c], h));
        }
      }
      using T = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<T>(&dl_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string fn = GetDumpName("s", ".vtk", par->dmp->GetN());
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " t=" << owner_->GetTime() + owner_->GetTimeStep()
            << " to " << fn << std::endl;
        WriteVtkPoly(dl_, fn);
      }
    }
  }
  void DumpParticles() {
    auto sem = m.GetSem("partdump");
    auto fr = par->part_dump_fr;
    size_t it = 1;
    if (1) { // TODO: revise frames
      if (sem("local")) {
        dpx_.clear();
        dpc_.clear();
        dpk_.clear();

        // copy to arrays  
        auto& bc = m.GetIndexCells();
        for (size_t s = 0; s < partstr_->GetNumStr(); ++s) {
          // cell
          IdxCell c = vsc_[s];
          auto v = GetPlaneBasis(m.GetCenter(c), fcn_[c], fca_[c], vsan_[s]);

          auto p = partstr_->GetStr(s);
          Vect* xx = p.first;
          size_t sx = p.second;

          auto w = bc.GetMIdx(c);
          const size_t mn = 1000; 
          // XXX: adhoc, hash for cell index, assume mesh size <= mn
          size_t ic = (w[2] * mn + w[1]) * mn + w[0]; // cell index
          for (size_t i = 0; i < sx; ++i) {
            auto x = GetSpaceCoords(xx[i], v);
            dpx_.push_back(x);
            dpc_.push_back(ic);
            dpk_.push_back(partstr_->GetCurv(s));
          }
        }

        // comm
        using TV = typename M::template OpCatT<Vect>;
        using TI = typename M::template OpCatT<size_t>;
        using TS = typename M::template OpCatT<Scal>;
        m.Reduce(std::make_shared<TV>(&dpx_));
        m.Reduce(std::make_shared<TI>(&dpc_));
        m.Reduce(std::make_shared<TS>(&dpk_));
      }
      if (sem("write")) {
        if (m.IsRoot()) {
          std::string s = GetDumpName("partit", ".csv", par->dmp->GetN(), 
                                      fr > 1 ? it : -1);
          std::cout << std::fixed << std::setprecision(8)
              << "dump" 
              << " t=" << owner_->GetTime() + owner_->GetTimeStep()
              << " to " << s << std::endl;
          std::ofstream o;
          o.open(s);
          o.precision(20);
          o << "x,y,z,c,k\n";

          for (size_t i = 0; i < dpx_.size(); ++i) {
            Vect x = dpx_[i];
            o << x[0] << "," << x[1] << "," << x[2] 
                << "," << dpc_[i] 
                << "," << dpk_[i] 
                << "\n";
          }
        }
      }
    }
  }
  Vect GetCellSize() const {
    Vect h; // cell size
    // XXX: specific for structured 3D mesh
    IdxCell c0(0);
    h = m.GetNode(m.GetNeighbourNode(c0, 7)) - 
        m.GetNode(m.GetNeighbourNode(c0, 0));
    assert(std::abs(h.prod() - m.GetVolume(c0)) < 1e-10);
    return h;
  }
  // Local coordinates in plane section of the interface.
  // xc: cell center
  // n: normal to interface
  // a: plane constant
  // an: angle
  // Returns:
  // {mc,mx,my,mn}: 
  // mc: coordinate center
  // mx: unit in x, tangent to interface
  // my: unit in y, normal to interface
  // mn: unit vector to from orthonormal positively oriented basis <mx,my,mn> 
  std::array<Vect, 3> GetPlaneBasis(const Vect& xc, const Vect& n, 
                                    Scal a, Scal an) {
    Vect h = GetCellSize();
    // direction in which normal has minimal component
    size_t d = n.abs().argmin(); 
    Vect xd(0); 
    xd[d] = 1.;
    // t0 orthogonal to n and d, <n,d,t0> positively oriented
    Vect t0 = n.cross(xd); 
    t0 /= t0.norm();
    // t1 orthogonal to n and t0, <t0,t1,n> positively oriented
    Vect t1 = n.cross(t0);
    t1 /= t1.norm();
    // tangent at angle an
    Vect t = t0 * std::cos(an) + t1 * std::sin(an); 
    t /= t.norm();

    Vect mc = xc + R::GetCenter(n, a, h); // center of interface
    Vect mx = t;  // unit in x, tangent
    Vect my = n / n.norm();  // unit in y, normal 
    return {mc, mx, my};
  }
  // Convert from space to plane coordinates.
  // x: space coordinates
  // v: output of GetPlaneBasis()
  // Returns:
  // Vect(xl,yl,0) with plane coordinates xl,yl
  Vect GetPlaneCoords(const Vect& x, const std::array<Vect, 3>& v) {
    auto mc = v[0];
    auto mx = v[1];
    auto my = v[2];
    return Vect((x - mc).dot(mx), (x - mc).dot(my), 0.);
  }
  // Convert from plane to space coordinates.
  // (xl,yl,0): plane coordinates
  // v: output of GetPlaneBasis()
  // Returns:
  // x: space coordinates
  Vect GetSpaceCoords(const Vect& xl, const std::array<Vect, 3>& v) {
    auto mc = v[0];
    auto mx = v[1];
    auto my = v[2];
    return mc + mx * xl[0] + my * xl[1];
  }
  void Seed0(const FieldCell<Scal>& fcu, const FieldCell<Scal>& fca, 
             const FieldCell<Vect>& fcn, const FieldCell<bool>& fci) {
    using MIdx = typename M::MIdx;
    auto& bc = m.GetIndexCells();
    Vect h = GetCellSize(); // cell size

    // clear string list
    partstr_->Clear();
    vsc_.clear();
    vsan_.clear();

    std::vector<std::array<Vect, 2>> ll; // buffer for interface lines
    // Seed strings in cells with interface.
    for (auto c : m.Cells()) {
      Vect xc = m.GetCenter(c);
      const Scal th = par->part_intth;
      if (fci[c] && fcu[c] >= th && fcu[c] <= 1. - th) {
        // number of strings
        size_t ns = (par->dim == 2 ? 1 : par->part_ns);
        for (size_t s = 0; s < ns; ++s) {
          Scal an = s * M_PI / ns; // angle
          auto v = GetPlaneBasis(xc, fcn[c], fca[c], an);

          // block of offsets to neighbours in stencil [-sw,sw]
          const int sw = 2; // stencil halfwidth
          const int sn = sw * 2 + 1; // stencil size
          GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, par->dim == 2 ? 0 : -sw), 
                                  MIdx(sn, sn, par->dim == 2 ? 1 : sn)); 

          auto w = bc.GetMIdx(c);
          ll.clear();
          // Extract interface from neighbour cells.
          for (auto wo : bo) {
            IdxCell cc = bc.GetIdx(w + wo); // neighbour cell
            if (fci[cc] && fcu[cc] >= th && fcu[cc] <= 1. - th) {
              auto xcc = m.GetCenter(cc); // center of cell
              auto xx = R::GetCutPoly(xcc, fcn[cc], fca[cc], h); // polygon
              std::array<Vect, 2> e; // ends of intersection 
              Vect mc = v[0];  // plane center
              Vect mx = v[1];  // unit in x
              Vect my = v[2];  // unit in y
              Vect mn = mx.cross(my); // normal to plane
              if (R::GetInterPoly(xx, mc, mn, e)) { // non-empty
                // interface normal 
                auto pncc = GetPlaneCoords(mc + fcn[cc], v);
                // line ends
                auto pe0 = GetPlaneCoords(e[0], v);
                auto pe1 = GetPlaneCoords(e[1], v);
                // make <pncc,pe1-pe0> positively oriented
                if (pncc.cross_third(pe1 - pe0) < 0.) {
                  std::swap(pe0, pe1);
                }
                ll.push_back({pe0, pe1});
              }
            }
          }

          // add string 
          partstr_->Add(Vect(0.), Vect(1., 0., 0.), ll.data(), ll.size());
          vsc_.push_back(c);
          vsan_.push_back(an);
          assert(vsc_.size() == partstr_->GetNumStr());
          assert(vsan_.size() == partstr_->GetNumStr());
        }
      }
    }
  }
  std::pair<Vect, Vect> GetBubble() {
    static Vect c(0);
    static Vect r(0);
    static bool ld = false;
    if (!ld) {
      std::ifstream f("../b.dat");
      f >> c[0] >> c[1] >> c[2];
      f >> r[0];
      f >> r[1];
      if (!f.good()) {
        r[1] = r[0];
      } else {
        f >> r[2];
        if (!f.good()) {
          r[2] = r[0];
        }
      }
      std::cout << "Loaded c=" << c << " r=" << r << std::endl;
      ld = true;
    }
    return std::make_pair(c, r);
  }
  void Part(const FieldCell<Scal>& uc, typename M::Sem& sem) {
    if (sem("part-comm")) {
      m.Comm(&fca_);
      m.Comm(&fcn_);
    }

    bool dm = par->dmp->Try(owner_->GetTime() + owner_->GetTimeStep(), 
                            owner_->GetTimeStep());

    if (sem("part-run")) {
      Seed0(uc, fca_, fcn_, fci_);
      partstr_->Run(par->part_tol, par->part_itermax, 
                    par->part_verb && m.IsRoot());
      // compute curvature
      fckp_.Reinit(m, GetNan<Scal>());
      // XXX: assume strings from same cell contiguous
      auto ns = partstr_->GetNumStr();
      size_t s = 0;
      while (s < ns) { // new cell
        size_t nsc = 0; // number of strings in current cell
        IdxCell c = vsc_[s];
        fckp_[c] = 0.; // reset sum
        while (s < ns && c == vsc_[s]) {
          fckp_[c] += partstr_->GetCurv(s);
          ++nsc;
          ++s;
        }
        fckp_[c] /= nsc;
        if (par->dim == 3) {
          fckp_[c] *= 2.;
        }
      }
      m.Comm(&fckp_);
    }
    if (sem("part-nan")) {
      // if interface cell but still nan, find from neighbour
      // block of offsets to neighbours in stencil [-sw,sw]
      const int sw = 1; // stencil halfwidth
      const int sn = sw * 2 + 1; // stencil size
      using MIdx = typename M::MIdx;
      GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, par->dim == 2 ? 0 : -sw), 
                              MIdx(sn, sn, par->dim == 2 ? 1 : sn)); 
      auto& bc = m.GetIndexCells();
      size_t nan = 0;
      IdxCell cnan;
      for (auto c : m.Cells()) {
        if (fci_[c] && IsNan(fckp_[c])) {
          auto w = bc.GetMIdx(c);
          for (auto wo : bo) {
            auto cc = bc.GetIdx(w + wo);
            if (!IsNan(fckp_[cc])) {
              fckp_[c] = fckp_[cc];
              break;
            }
          }
          if (IsNan(fckp_[c])) { // still nan, fill with zero
            ++nan;
            cnan = c;
            fckp_[c] = 0.;
          }
        }
      }
      if (nan && par->part_verb) {
        std::stringstream s;
        s.precision(16);
        s << "nan curvature in all neighbours of " << nan 
            << " cells, one at x=" << m.GetCenter(cnan)
            << " with vf=" << uc[cnan];
        //throw std::runtime_error(s.str());
        std::cout << s.str() << std::endl;
      }
      m.Comm(&fckp_);
    }

    if (sem.Nested("part-dump")) {
      if (dm) {
        DumpParticles();
      }
    }
  }

  // reconstruct interface
  void Rec(const FieldCell<Scal>& uc) {
    auto sem = m.GetSem("reconst");

    // XXX: adhoc
    // overwrite u=0 if y<y0 or y>y0
    if (sem("bc-zero")) {
      Scal y0 = par->bcc_y0;
      Scal y1 = par->bcc_y1;
      auto& uu = const_cast<FieldCell<Scal>&>(uc);
      for (auto c : m.AllCells()) {
        auto x = m.GetCenter(c);
        if (x[1] < y0 || x[1] > y1) {
          uu[c] = 0.;
        }
      }
    }

    if (sem("height")) {
      CheckNan(uc, "vof:Rec:uc", m);
      DetectInterface(uc);
      // Compute normal and curvature [s]
      CalcNormal(uc, fci_, fcn_, fck_);
      auto h = GetCellSize();
      // Reconstruct interface [s]
      for (auto c : m.SuCells()) {
        fca_[c] = R::GetLineA(fcn_[c], uc[c], h);
      }
    }

    // Correction with normal from particles
    if (par->part && par->part_n) {
      Part(uc, sem);
      if (sem("parta")) {
        auto h = GetCellSize();
        for (auto c : m.AllCells()) {
          fca_[c] = R::GetLineA(fcn_[c], uc[c], h);
        }
      }
    }
  }
  void CalcNormal(const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
                  FieldCell<Vect>& fcn, FieldCell<Scal>& fck) {
    UNormal<M>::CalcNormal(m, fcu, fci, par->dim, fcn, fck);
  }
  void DetectInterface(const FieldCell<Scal>& uc) {
    fci_.Reinit(m, false);
    // volume fraction different from 0 or 1
    for (auto c : m.AllCells()) {
      Scal u = uc[c];
      if (u > 0. && u < 1.) {
        fci_[c] = true;
      }
    }
    // neighbour cell has different value but both are 0 or 1
    for (auto f : m.SuFaces()) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      Scal um = uc[cm];
      Scal up = uc[cp];
      if ((um == 0. || um == 1.) && (up == 0. || up == 1.) && (um != up)) {
        fci_[cm] = true;
        fci_[cp] = true;
      }
    }
  }
  void StartStep() {
    auto sem = m.GetSem("start");
    if (sem("rotate")) {
      owner_->ClearIter();
      fcu_.time_prev = fcu_.time_curr;
      fcu_.iter_curr = fcu_.time_prev;
    }

    if (sem.Nested("reconst")) {
      if (owner_->GetTime() == 0.) {
        Rec(fcu_.time_curr);
      }
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      auto& uc = fcu_.iter_curr;
      const Scal dt = owner_->GetTimeStep();
      auto& fcs = *owner_->fcs_;
      for (auto c : m.Cells()) {
        uc[c] = 
            fcu_.time_prev[c] +  // previous time step
            dt * fcs[c]; // source
      }
    }

    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    // directions, format: {dir LE, dir EI, ...}
    std::vector<size_t> dd; 
    Scal vsc; // scaling factor for ffv, used for splitting
    if (par->dim == 3) { // 3d
      if (count_ % 3 == 0) {
        dd = {0, 1, 1, 2, 2, 0};
      } else if (count_ % 3 == 1) {
        dd = {1, 2, 2, 0, 0, 1};
      } else {
        dd = {2, 0, 0, 1, 1, 2};
      }
      vsc = 0.5;
    } else {
      if (count_ % 2 == 0) {
        dd = {0, 1};
      } else {
        dd = {1, 0};
      } 
      vsc = 1.0;
    }
    for (size_t id = 0; id < dd.size(); ++id) {
      if (sem("adv")) {
        size_t d = dd[id]; // direction as index
        Dir md(d); // direction as Dir
        MIdx wd(md); // offset in direction d
        auto& uc = fcu_.iter_curr;
        auto& bc = m.GetIndexCells();
        auto& bf = m.GetIndexFaces();
        auto h = GetCellSize();
        auto& ffv = *owner_->ffv_; // [f]ield [f]ace [v]olume flux
        const Scal dt = owner_->GetTimeStep();

        FieldFace<Scal> ffvu(m); // flux: volume flux * field

        for (auto f : m.Faces()) {
          auto p = bf.GetMIdxDir(f);
          MIdx wf = p.first;
          Dir df = p.second;

          if (df != md) {
            continue;
          }

          // mixture flux
          const Scal v = ffv[f] * vsc;
          // upwind cell
          IdxCell cu = m.GetNeighbourCell(f, v > 0. ? 0 : 1);
          if (fci_[cu]) { // cell contains interface, flux from reconstruction
            if (id % 2 == 0) { // Euler Implicit
              // phase 2 flux
              ffvu[f] = R::GetLineFlux(fcn_[cu], fca_[cu], h, v, dt, d);
            } else { // Lagrange Explicit
              // upwind face
              IdxFace fu = bf.GetIdx(v > 0. ? wf - wd : wf + wd, md);
              // upwind mixture flux
              Scal vu = ffv[fu] * vsc;
              // phase 2 flux
              ffvu[f] = R::GetLineFluxStr(fcn_[cu], fca_[cu], h, v, vu, dt, d);
            }
          } else {
            ffvu[f] = v * uc[cu];
          }
        }

        FieldFace<Scal> ffu(m);
        // interpolate field value to boundaries
        InterpolateB(uc, mfc_, ffu, m);

        // override boundary upwind flux
        for (const auto& it : mfc_) {
          IdxFace f = it.GetIdx();
          CondFace* cb = it.GetValue().get(); 
          Scal v = ffv[f];
          if ((cb->GetNci() == 0) != (v > 0.)) {
            ffvu[f] = v * ffu[f];
            // XXX: adhoc
            // Alternating mul correction of flux
            // (done for bubble detachment)
            if (par->bcc_t0 > 0. && par->bcc_t1 > 0.) {
              Scal k0 = par->bcc_k0;
              Scal k1 = par->bcc_k1;
              Scal t0 = par->bcc_t0;
              Scal t1 = par->bcc_t1;
              Scal t = owner_->GetTime();
              Scal ts = t0 + t1;
              Scal ph = t / ts; 
              ph = ph - int(ph);
              ph *= ts;
              ffvu[f] *= (ph < t0 ? k0 : k1);
            }
          }
        }

        for (auto c : m.Cells()) {
          auto w = bc.GetMIdx(c);
          const Scal lc = m.GetVolume(c);
          // faces
          IdxFace fm = bf.GetIdx(w, md);
          IdxFace fp = bf.GetIdx(w + wd, md);
          // mixture fluxes
          const Scal vm = ffv[fm] * vsc;
          const Scal vp = ffv[fp] * vsc;
          // mixture cfl
          const Scal sm = vm * dt / lc;
          const Scal sp = vp * dt / lc;
          const Scal ds = sp - sm;
          // phase 2 fluxes
          Scal qm = ffvu[fm];
          Scal qp = ffvu[fp];
          // phase 2 cfl
          const Scal lm = qm * dt / lc;
          const Scal lp = qp * dt / lc;
          const Scal dl = lp - lm;
          if (id % 2 == 0) { // Euler Implicit
            uc[c] = (uc[c] - dl) / (1. - ds);
          } else { // Lagrange Explicit
            uc[c] = uc[c] * (1. + ds) - dl;
          }
        }

        // clip
        const Scal th = par->clipth;
        for (auto c : m.Cells()) {
          Scal& u = uc[c];
          if (u < th) {
            u = 0.;
          } else if (u > 1. - th) {
            u = 1.;
          } else if (IsNan(u)) {
            u = 0.;
          }
        }
        m.Comm(&uc);
      }
      if (sem.Nested("reconst")) {
        Rec(fcu_.iter_curr);
      }
    }

    // Curvature from gradient of volume fraction
    if (par->curvgrad && sem("curv")) {
      auto ffu = Interpolate(fcu_.iter_curr, mfc_, m); // [s]
      auto fcg = Gradient(ffu, m); // [s]
      auto ffg = Interpolate(fcg, mfvz_, m); // [i]

      fck_.Reinit(m, GetNan<Scal>()); // curvature [i]
      for (auto c : m.Cells()) {
        if (!fci_[c]) {
          continue;
        }
        Scal s = 0.;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          auto& g = ffg[f];
          auto n = g / g.norm();  // inner normal
          s += -n.dot(m.GetOutwardSurface(c, q));
        }
        fck_[c] = s / m.GetVolume(c);
      }
    }

    if (sem("curvcomm")) {
      m.Comm(&fck_);
    }

    if (par->dumppoly) {
      bool dm = par->dmp->Try(owner_->GetTime() + owner_->GetTimeStep(), 
                              owner_->GetTimeStep());
      if (dm && sem("dumppoly")) {
        DumpPoly();
      }
    }

    if (par->part) {
      Part(fcu_.iter_curr, sem);
    }

    if (sem("stat")) {
      owner_->IncIter();
      ++count_;
    }
  }
  void FinishStep() {
    fcu_.time_curr = fcu_.iter_curr;
    owner_->IncTime();
  }

  Vof* owner_;
  std::shared_ptr<Par> par;
  M& m;

  LayersData<FieldCell<Scal>> fcu_;
  MapFace<std::shared_ptr<CondFace>> mfc_;
  MapFace<std::shared_ptr<CondFace>> mfvz_; // zero-derivative bc for Vect

  FieldCell<Scal> fca_; // alpha (plane constant)
  FieldCell<Vect> fcn_; // n (normal to plane)
  FieldCell<Scal> fck_; // curvature from height functions
  FieldCell<Scal> fckp_; // curvature from particles
  FieldCell<bool> fci_; // interface mask (1: contains interface)
  size_t count_ = 0; // number of MakeIter() calls, used for splitting
  // XXX: fca_ and fcn_ are proportional, separate scaling is invalid

  std::unique_ptr<PS> partstr_; // particle strings
  std::vector<IdxCell> vsc_; // vsc_[s] is cell of string s
  std::vector<Scal> vsan_; // vsc_[s] is angle of tangent (see GetPlaneBasis)

  std::vector<Vect> dpx_; // dump particles x
  std::vector<size_t> dpc_; // dump particles cell
  std::vector<Scal> dpk_; // dump particles curvature
  std::vector<std::vector<Vect>> dl_; // dump poly
};

template <class M_>
Vof<M_>::Vof(M& m, const FieldCell<Scal>& fcu,
    const MapFace<std::shared_ptr<CondFace>>& mfc,
    const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
    double t, double dt, std::shared_ptr<Par> par)
    : AdvectionSolver<M>(t, dt, m, ffv, fcs)
    , imp(new Imp(this, fcu, mfc, par))
{}

template <class M_>
auto Vof<M_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_>
void Vof<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void Vof<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void Vof<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
auto Vof<M_>::GetField(Layers l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}

template <class M_>
auto Vof<M_>::GetAlpha() const -> const FieldCell<Scal>& {
  return imp->fca_;
}

template <class M_>
auto Vof<M_>::GetNormal() const -> const FieldCell<Vect>& {
  return imp->fcn_;
}

template <class M_>
auto Vof<M_>::GetCurv() const -> const FieldCell<Scal>& {
  return imp->par->part_k ? imp->fckp_ : imp->fck_;
}

// curvature from height function
template <class M_>
auto Vof<M_>::GetCurvH() const -> const FieldCell<Scal>& {
  return imp->fck_;
}

// curvature from particles
template <class M_>
auto Vof<M_>::GetCurvP() const -> const FieldCell<Scal>& {
  return imp->fckp_;
}

} // namespace solver
