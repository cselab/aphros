#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>

#include "advection.h"
#include "geom/block.h"
#include "dump/dumper.h"
#include "dump/vtk.h"
#include "reconst.h"
#include "partstr.h"

// attraction to exact sphere
#define ADHOC_ATTR 0
// normal from exact sphere
#define ADHOC_NORM 0

namespace solver {

template <class M_>
class Vof : public AdvectionSolver<M_> {
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using R = Reconst<Scal>;
  using PS = PartStr<Scal>;
  static constexpr size_t dim = M::dim;

  using P::m;
  using P::ffv_;
  using P::fcs_;
  LayersData<FieldCell<Scal>> fc_u_;
  MapFace<std::shared_ptr<CondFace>> mfc_;
  MapFace<std::shared_ptr<CondFace>> mfvz_; // zero-derivative bc for Vect

  FieldCell<Scal> fc_a_; // alpha (plane constant)
  FieldCell<Vect> fc_n_; // n (normal to plane)
  FieldCell<Scal> fc_us_; // smooth field
  FieldFace<Scal> ff_fu_; // volume flux
  FieldCell<Scal> fck_; // curvature from height functions
  FieldCell<Scal> fckp_; // curvature from particles
  FieldFace<Scal> ffu_; // field on faces
  FieldFace<Scal> ffvu_; // flux: volume flux * field
  FieldCell<Vect> fcg_; // gradient
  FieldFace<Vect> ffg_; // gradient
  FieldCell<bool> fci_; // interface mask (1: contains interface)
  size_t count_ = 0; // number of MakeIter() calls, used for splitting
  static constexpr size_t kNp = 11; // particles in one string
  static constexpr size_t kNpp = 4; // maximum strings per cell
  FieldCell<std::array<Vect, kNp * kNpp>> fcp_; // cell list
  FieldCell<std::array<Vect, kNp * kNpp>> fcpt_; // cell list tmp
  FieldCell<std::array<Scal, kNp * kNpp>> fcpw_; // cell list weight
  FieldCell<size_t> fcps_; // cell list size

  using PSP = typename PS::Par;
  std::unique_ptr<PS> partstr_; // particle strings
  std::vector<IdxCell> vsc_; // vsc_[s] is cell of string s

  std::vector<Vect> dpx_; // dump particles x
  std::vector<size_t> dpc_; // dump particles cell
  std::vector<std::vector<Vect>> dl_; // dump poly

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
  }

 public:
  struct Par {
    size_t dim = 3; // dimension (dim=2 assumes zero velocity in z)
    bool curvgrad = false; // compute curvature using gradient
    bool part = false; // particles
    Scal part_relax = 1.; 
    Scal part_h0 = 1.; // dist init
    Scal part_h = 1.;  // dist eq
    Scal part_kstr = 1.; // stretching
    Scal part_kattr = 1.; // attraction to reconstructed interface
    Scal part_kbend = 1.; // bending
    bool part_bendmean = true; // bending to mean angle (fit circle)
    bool part_n = false; // normal from particles
    // curvature from particles
    // if true, GetCurv returns fckp_
    bool part_k = false; // curvature from particles
    size_t part_maxiter = 100; // num iter
    size_t part_dump_fr = 100; // num frames dump
    size_t part_report_fr = 100; // num frames report
    Scal part_intth = 1e-3; // interface detection threshold
    Scal clipth = 1e-6; // vf clipping threshold
    std::unique_ptr<Dumper> dmp; // dumper for particles
    bool dumppoly = false; // dump reconstructed interface (cut polygons)
    // XXX: adhoc
    Scal bcc_k0 = 1.;   // mul corrections to bc 
    Scal bcc_k1 = 1.;   
    Scal bcc_t0 = -1.;   // duration of phases (one negative to disable)
    Scal bcc_t1 = -1.;   
    Scal bcc_y0 = -1e10; // overwrite u=0 if y<y0 or y>y1
    Scal bcc_y1 = 1e10;  // (to remove periodic conditions)
    int part_constr = 0; // 0: no constraints
                         // 1: fixed distance, constant angle
                         // 2: fixed distance, linear angle
    Scal part_segcirc = 1.; // factor for circular segment
    size_t part_np = 11; // number of particles per string
    size_t part_ns = 4; // number of strings per cell
    size_t part_itermax = 100; // itermax
    Scal part_tol = 0.01; // tolerance
  };
  std::shared_ptr<Par> par;
  Par* GetPar() { return par.get(); }
  void SeedParticles(const FieldCell<Scal>& uc) {
    fcps_.Reinit(m, 0);
    fcp_.Reinit(m);
    fcpt_.Reinit(m);
    fcpw_.Reinit(m);

    auto& bc = m.GetBlockCells();
    auto& bn = m.GetBlockNodes();
    using MIdx = typename M::MIdx;
    MIdx wb = bn.GetBegin();
    Vect xb = m.GetNode(bn.GetIdx(wb));
    Vect h = m.GetNode(bn.GetIdx(wb + MIdx(1))) - xb;
    Scal hm = h.norminf();

    for (auto c : m.Cells()) {
      const Scal th = par->part_intth;
      if (fci_[c] && uc[c] > th && uc[c] < 1. - th) {
        Vect x = m.GetCenter(c) + R::GetCenter(fc_n_[c], fc_a_[c], h);
        Vect n = fc_n_[c];
        n /= n.norm();
        // direction in which normal has minimal component
        size_t d = n.abs().argmin(); 
        Vect xd(0); 
        xd[d] = 1.;
        // t0 orthogonal to n and d, <n,d,t0> positively oriented
        Vect t0 = n.cross(xd); 
        t0 /= t0.norm();
        // t1 orthogonal to n and t0
        Vect t1 = n.cross(t0);
        t1 /= t1.norm();
        const Scal pd = hm * par->part_h0; // distance between particles
        if (fcps_[c] == 0) { // if no particles yet
          Scal a = 0.; // angle of tangent
          // number of strings
          size_t ns = (par->dim == 2 ? 1 : kNpp);
          for (int s = 0; s < ns; ++s) {
            // tangent to string
            Vect t = t0 * std::cos(a) + t1 * std::sin(a);
            for (int i = 0; i < kNp; ++i) {
              fcp_[c][fcps_[c]++] = x + t * ((i - (kNp - 1) * 0.5) * pd);
            }
            if (ns > 1) {
              a += M_PI / ns;
            }
          }
        }
      }
    }
  }
  // Dump cut polygons
  void DumpPoly() {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      auto h = GetCellSize();
      for (auto c : m.Cells()) {
        const Scal th = par->part_intth;
        Scal u = fc_u_.iter_curr[c];
        if (fci_[c] && u > th && u < 1. - th) {
          dl_.push_back(R::GetCutPoly(m.GetCenter(c), fc_n_[c], fc_a_[c], h));
        }
      }
      using T = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<T>(&dl_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string st = "." + std::to_string(par->dmp->GetN());
        auto fn = "s" + st + ".vtk";
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " t=" << this->GetTime() + this->GetTimeStep()
            << " to " << fn << std::endl;
        WriteVtkPoly(dl_, fn);
      }
    }
  }
  void DumpParticles(size_t it) {
    auto sem = m.GetSem("partdump");

    auto fr = par->part_dump_fr;
    size_t d = std::max<size_t>(1, par->part_maxiter / fr);
    if (fr > 1 && it % d == 0 || it + 1 == par->part_maxiter) {
      if (sem("local")) {
        dpx_.clear();
        dpc_.clear();

        // copy to arrays  
        auto& bc = m.GetBlockCells();
        for (auto c : m.Cells()) {
          for (size_t i = 0; i < fcps_[c]; ++i) {
            dpx_.push_back(fcp_[c][i]);
            auto w = bc.GetMIdx(c);
            // XXX: adhoc, hash for cell index, assume mesh size <= mn
            const size_t mn = 1000; 
            dpc_.push_back((w[2] * mn + w[1]) * mn + w[0]);
          }
        }

        // comm
        using TV = typename M::template OpCatT<Vect>;
        using TI = typename M::template OpCatT<size_t>;
        m.Reduce(std::make_shared<TV>(&dpx_));
        m.Reduce(std::make_shared<TI>(&dpc_));
      }
      if (sem("write")) {
        if (m.IsRoot()) {
          std::string st = "." + std::to_string(par->dmp->GetN());
          std::string sit = fr > 1 ? "_" + std::to_string(it) : "";
          std::string s = "partit" + st + sit + ".csv";
          std::cout << std::fixed << std::setprecision(8)
              << "dump" 
              << " t=" << this->GetTime() + this->GetTimeStep()
              << " to " << s << std::endl;
          std::ofstream o;
          o.open(s);
          o.precision(20);
          o << "x,y,z,c\n";

          for (size_t i = 0; i < dpx_.size(); ++i) {
            Vect x = dpx_[i];
            o << x[0] << "," << x[1] << "," << x[2] 
                << "," << dpc_[i] << "\n";
          }
        }
      }
    }
  }
  Vof(M& m, const FieldCell<Scal>& fcu,
      const MapFace<std::shared_ptr<CondFace>>& mfc,
      const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
      double t, double dt, std::shared_ptr<Par> par)
      : AdvectionSolver<M>(t, dt, m, ffv, fcs)
      , mfc_(mfc), par(par)
      , fc_a_(m, 0), fc_n_(m, Vect(0)), fc_us_(m, 0), ff_fu_(m, 0) 
      , fck_(m, 0), fckp_(m, 0)
  {
    fc_u_.time_curr = fcu;
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
  void StartStep() override {
    auto sem = m.GetSem("start");
    if (sem("rotate")) {
      this->ClearIter();
      fc_u_.time_prev = fc_u_.time_curr;
      fc_u_.iter_curr = fc_u_.time_prev;
    }

    if (sem.Nested("reconst")) {
      if (this->GetTime() == 0.) {
        Rec(fc_u_.time_curr);
      }
    }

      /*
      if (par->part) {
        if (sem("seed")) {
          SeedParticles(fc_u_.time_curr);
          if (par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                            this->GetTimeStep())) {
            DumpParticles(par->part_maxiter - 1);
          }
        }
      }
      */
  }
  Scal Maxmod(Scal a, Scal b) {
    return std::abs(b) < std::abs(a) ? a : b;
  }
  // Computes normal by gradient.
  // uc: volume fraction
  // msk: interface mask (1: contains interface)
  // Output: modified in cells with msk=1
  // fcn: normal [s] 
  void CalcNormalGrad(const FieldCell<Scal>& uc) {
    auto uf = Interpolate(uc, mfc_, m);
    auto gc = Gradient(uf, m);
    for (auto c : m.AllCells()) {
      Vect g = gc[c];
      fc_n_[c] = g;
    }
  }
  // Computes normal by Young's scheme (interpolation from nodes).
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // Output: modified in cells with fci=1, resized to m
  // fcn: normal with norm1()=1, antigradient of fcu [s] 
  // XXX: uses static variables, not suspendable
  // TODO: check non-uniform mesh
  void CalcNormalYoung(const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
                       FieldCell<Vect>& fcn) {
    static FieldNode<Vect> g;  // gradient
    static FieldNode<Vect> l;  // step from cell to node
    g.Reinit(m, Vect(0));
    l.Reinit(m, Vect(0));
    // values from cells to neighbour nodes
    for (auto c : m.AllCells()) {
      Vect xc = m.GetCenter(c);
      for (size_t q = 0; q < m.GetNumNeighbourNodes(c); ++q) {
        IdxNode n = m.GetNeighbourNode(c, q);
        Vect xn = m.GetNode(n);
        for (size_t d = 0; d < dim; ++d) {
          g[n][d] += (xc[d] - xn[d] > 0. ? 1. : -1.) * fcu[c];
          l[n][d] += std::abs(xc[d] - xn[d]);
        }
      } 
    } 
    // gradient on nodes
    for (auto n : m.SuNodes()) {
      g[n] /= l[n];
    }

    // gradient on cells
    fcn.Reinit(m);
    for (auto c : m.SuCells()) {
      if (fci[c]) {
        // sum over neighbour nodes
        auto& v = fcn[c];
        v = Vect(0);
        for (size_t q = 0; q < m.GetNumNeighbourNodes(c); ++q) {
          IdxNode n = m.GetNeighbourNode(c, q);
          v += g[n];
        }
        // normalize
        v /= -v.norm1();
      }
    }
  }
  // Computes normal and curvature from height functions.
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // ow: 1: force overwrite, 0: update only if gives steeper profile
  // Output: modified in cells with fci=1, resized to m
  // fcn: normal, antigradient of fcu, if gives steeper profile or ow=1 [s] 
  // fck: curvature [s] 
  void CalcNormalHeight(const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
                        bool ow, FieldCell<Vect>& fcn, FieldCell<Scal>& fck) {
    using MIdx = typename M::MIdx;
    using Dir = typename M::Dir;
    auto& bc = m.GetBlockCells();

    fcn.Reinit(m); // XXX
    fck.Reinit(m);

    for (auto c : m.SuCells()) {
      if (!fci[c]) {
        continue;
      }
      Vect bn; // best normal
      Scal bhx, bhy; // best first derivative
      Scal bk; // best curvature[k]
      Dir bd;  // best direction
      bool fst = true; // first
      std::vector<Dir> dd; // direction of plane normal
      if (par->dim == 2) {
        dd = {Dir::i, Dir::j};
      } else {
        dd = {Dir::i, Dir::j, Dir::k};
      }
      for (Dir dn : dd) {
        // directions of plane tangents ([d]irection [t]angents)
        Dir dtx((size_t(dn) + 1) % dim); 
        Dir dty((size_t(dn) + 2) % dim); 

        MIdx w = bc.GetMIdx(c);

        // offset in normal direction
        MIdx on = MIdx(dn);
        // offset in dtx,dty
        MIdx otx = MIdx(dtx);
        MIdx oty = MIdx(dty);
        // mesh step
        const Scal lx = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - otx)));
        const Scal ly = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - oty)));
        const Scal ln = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - on)));

        // Evaluates height function
        // o: offset from w
        auto hh = [&](MIdx o) -> Scal {
          return 
            (fcu[bc.GetIdx(w + o - on)] + 
            fcu[bc.GetIdx(w + o)] + 
            fcu[bc.GetIdx(w + o + on)]) * ln;
        };

        // height function
        const Scal h = hh(MIdx(0));
        const Scal hxm = hh(-otx);
        const Scal hxp = hh(otx);
        const Scal hym = hh(-oty);
        const Scal hyp = hh(oty);
        // corners: hxy
        const Scal hmm = hh(-otx - oty); 
        const Scal hmp = hh(-otx + oty);
        const Scal hpm = hh(otx - oty);
        const Scal hpp = hh(otx + oty);

        // first derivative (slope)
        Scal hx = (hxp - hxm) / (2. * lx);  // centered
        Scal hy = (hyp - hym) / (2. * ly); 
        // sign: +1 if u increases in dn
        Scal sg = 
            (fcu[bc.GetIdx(w + on)] - fcu[bc.GetIdx(w - on)] > 0. ? 1. : -1.);
        // second derivative 
        Scal hxx = (hxp - 2. * h + hxm) / (lx * lx);
        Scal hyy = (hyp - 2. * h + hym) / (ly * ly);
        Scal hxy = ((hpp - hmp) - (hpm - hmm)) / (4. * lx * ly);
        // curvature
        Scal k = (2. * hx * hy * hxy 
            -(sqr(hy) + 1.) * hxx -(sqr(hx) + 1.) * hyy) / 
            std::pow(sqr(hx) + sqr(hy) + 1., 3. / 2.);
        // outer normal
        Vect n;
        n[size_t(dtx)] = -hx;
        n[size_t(dty)] = -hy;
        n[size_t(dn)] = -sg;
        // select best with minimal slope
        if (fst || 
            std::abs(hx) + std::abs(hy) < std::abs(bhx) + std::abs(bhy)) {
          bn = n;
          bhx = hx;
          bhy = hy;
          bk = k;
          bd = dn;
          fst = false;
        } 
      }
      bn /= bn.norm1(); // normalize

      // update if ow=1 or gives steeper profile in plane dn
      if (ow || std::abs(bn[size_t(bd)]) < std::abs(fcn[c][size_t(bd)])) {
        fcn[c] = bn;
      }

      // curvature
      fck[c] = bk;

      #if ADHOC_NORM
      auto cr = GetBubble();
      auto q = m.GetCenter(c) - cr.first;
      if (par->dim == 2) {
        q[2] = 0.;
      }
      fcn[c] = q / q.norm();
      #endif 
    }
  }
  // Computes normal by combined Young scheme and height-functions
  // and curvature from height-functions.
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // Output: set to NaN if fci=0
  // fcn: normal with norm1()=1, antigradient of fcu [s] 
  // fck: curvature
  void CalcNormal(const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
                  FieldCell<Vect>& fcn, FieldCell<Scal>& fck) {
    fcn.Reinit(m, Vect(GetNan<Scal>()));
    fck.Reinit(m, GetNan<Scal>());
    CalcNormalYoung(fcu, fci, fcn);
    CalcNormalHeight(fcu, fci, false, fcn, fck);
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
  void Seed0(const FieldCell<Scal>& fcu, const FieldCell<Scal>& fca, 
             const FieldCell<Vect>& fcn, const FieldCell<bool>& fci) {
    using MIdx = typename M::MIdx;
    auto& bc = m.GetBlockCells();
    Vect h = GetCellSize(); // cell size

    std::vector<std::array<Vect, 2>> ll; // buffer for interface lines
    for (auto c : m.Cells()) {
      const Scal th = par->part_intth;
      if (fci[c] && fcu[c] >= th && fcu[c] <= 1. - th) {
        // center of interface
        Vect xc = m.GetCenter(c) + R::GetCenter(fcn[c], fca[c], h);
        // unit normal to interface
        Vect n = fcn[c];
        n /= n.norm();
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

        // number of strings
        size_t ns = (par->dim == 2 ? 1 : par->part_ns);
        for (int s = 0; s < ns; ++s) {
          Scal a = s * M_PI / ns; // angle
          Vect t = t0 * std::cos(a) + t1 * std::sin(a); // tangent to string
          t /= t.norm();

          // Local coordinates in plane 
          Vect mx = t;  // unit in x
          Vect my = n;  // unit in y
          Vect mc = xc;  // center
          Vect mn = mx.cross(my);  // normal to plane
          mn /= mn.norm();

          // Transforms from space to plane coordinates.
          // x: space coordinates
          // Returns:
          // Vect(xl, yl, 0): plane coordinates
          auto pr = [&](Vect x) -> Vect {
            Scal xl = (x - mc).dot(mx);
            Scal yl = (x - mc).dot(my);
            return Vect(xl, yl, 0.);
          };

          // block of offsets to neighbours in stencil [-sw,sw]
          const int sw = 2; // stencil halfwidth
          const int sn = sw * 2 + 1; // stencil size
          GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, par->dim == 2 ? 0 : -sw), 
                                  MIdx(sn, sn, par->dim == 2 ? 1 : sn)); 

          auto w = bc.GetMIdx(c);
          // Cut interface in neighbour cells by plane.
          ll.clear();
          for (auto wo : bo) {
            IdxCell cc = bc.GetIdx(w + wo); // neighbour cell
            if (fci[cc] && fcu[cc] >= th && fcu[cc] <= 1. - th) {
              // center of cell
              auto xcc = m.GetCenter(cc);
              // interface polygon
              auto xx = R::GetCutPoly(xcc, fcn[cc], fca[cc], h);
              // ends of intersection 
              std::array<Vect, 2> e;
              if (R::GetInterPoly(xx, mc, mn, e)) { // non-empty
                // points in local coordinates
                // interface normal 
                auto pncc = pr(mc + fcn[cc]);
                // line ends
                auto pe0 = pr(e[0]);
                auto pe1 = pr(e[1]);
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
        }
      }
    }
  }
  // Compute force to advance particles.
  // Assume 2d positions (x,y,0).
  // x: pointer to first particle
  // sx: size of particle string
  // l: pointer to first array of line ends
  // nl: number of lines
  // Output:
  // f: position corrections of size sx
  void PartForce(const Vect* xx, size_t sx, 
                   const std::array<Vect, 2>* ll, size_t sl, Vect* ff) {
    PS::InterfaceForce(xx, sx, ll, sl, par->part_segcirc, ff);

    Scal hm = GetCellSize().norminf(); // cell size
    const int cn = par->part_constr; // constraint type
    if (cn == 0) {
      PS::Constr0(xx, sx, par->part_kstr, par->part_h * hm,
                  par->part_kbend, par->part_bendmean, par->part_relax, ff);
    } else if (cn == 1) {
      PS::Constr1(xx, sx, par->part_kattr, par->part_kbend, par->part_kstr,
                  hm, par->part_relax, ff);
    } else {
      throw std::runtime_error("Unknown part_constr=" + std::to_string(cn));
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
    if (sem("part-comma")) {
      m.Comm(&fc_a_);
      m.Comm(&fc_n_);
    }

    bool dm = par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                            this->GetTimeStep());

    if (sem("part-seed")) {
      SeedParticles(uc);

      Seed0(fc_u_.iter_curr, fc_a_, fc_n_, fci_);
      partstr_->Run(par->part_tol, par->part_itermax);
    }

    if (0 && sem.Nested("part-dump0")) {
      if (dm) {
        DumpParticles(0);
      }
    }

    if (0 && sem("part-advance")) {
      // particle strings
      std::vector<IdxCell> sc; // cell
      std::vector<Vect> xx; // particle positions
      std::vector<size_t> sx; // xx index, plus element sx.size() 
      std::vector<std::array<Vect, 2>> ll; // interface lines
      std::vector<size_t> sl; // sl index, plus element sl.size() 
      std::vector<Vect> smx; // local unit x
      std::vector<Vect> smy; // local unit y
      std::vector<Vect> smc; // local center
      std::vector<Scal> sk; // curvature

      // Extract interface, project particles
      for (auto c : m.Cells()) {
        // XXX: assume fcps_[c] % kNp == 0
        for (int is = 0; is < fcps_[c] / kNp; ++is) {
          int i0 = is * kNp;
          int i1 = (is + 1) * kNp;

          sc.push_back(c);

          // Plane coordinates
          // center
          Vect rc = fcp_[c][(i0 + i1 + 1) / 2];
          // tangent, assume straight line
          Vect rt = (fcp_[c][i0 + 1] - fcp_[c][i0]);
          rt /= rt.norm();
          // interface normal, assume orthogonal to rt
          Vect n = fc_n_[c];
          n /= n.norm();
          // string plane normal
          Vect rn = rt.cross(n);
          rn /= rn.norm();

          // unit in x
          Vect mx = rt;
          // unit in y
          Vect my = n;
          // center
          Vect mc = rc;

          smx.push_back(mx);
          smy.push_back(my);
          smc.push_back(mc);

          // Transform to plane coordinates.
          // x: space coordinates
          // Returns:
          // Vect(xl, yl, 0): plane coordinates
          auto pr = [&](Vect x) -> Vect {
            Vect q(0);
            q[0] = (x - mc).dot(mx);
            q[1] = (x - mc).dot(my);
            return q;
          };

          // Copy projected particles
          sx.push_back(xx.size());
          for (int i = i0; i < i1; ++i) {
            xx.push_back(pr(fcp_[c][i]));
          }

          // Extract interface lines
          sl.push_back(ll.size());
          auto& bc = m.GetBlockCells();
          using MIdx = typename M::MIdx;
          Vect h = GetCellSize();

          const int sw = 2; // stencil halfwidth, [-sw,sw]
          const int sn = sw * 2 + 1; // stencil size

          // block of offsets
          GBlock<IdxCell, dim> bo(
              MIdx(-sw, -sw, par->dim == 2 ? 0 : -sw), 
              MIdx(sn, sn, par->dim == 2 ? 1 : sn)); 

          auto w = bc.GetMIdx(c);
          for (auto wo : bo) {
            auto cc = bc.GetIdx(w + wo);
            Scal u = uc[cc];
            Scal th = par->part_intth;
            if (fci_[c] && u > th && u < 1. - th) {
              auto xcc = m.GetCenter(cc);
              auto xx = R::GetCutPoly(xcc, fc_n_[cc], fc_a_[cc], h);
              std::array<Vect, 2> e;
              if (R::GetInterPoly(xx, rc, rn, e)) {
                // projected line ends 
                // <pncc,pe1-pe0> positively oriented
                auto pncc = pr(mc + fc_n_[cc]);
                auto pe0 = pr(e[0]);
                auto pe1 = pr(e[1]);
                if (pncc.cross_third(pe1 - pe0) < 0.) {
                  std::swap(pe0, pe1);
                }
                ll.push_back({pe0, pe1});
              }
              #if ADHOC_ATTR 
              auto cr = GetBubble();
              ll.push_back({pr(cr.first), cr.second}); 
              #endif
            }
          }
        }
      }
      sx.push_back(xx.size());
      sl.push_back(ll.size());

      assert(sx.size() == sc.size() + 1);
      assert(sl.size() == sc.size() + 1);
      assert(smx.size() == sc.size());
      assert(smy.size() == sc.size());
      assert(smc.size() == sc.size());

      sk.resize(sc.size(), 0.);

      std::vector<Vect> ff(xx.size()); // force

      // advance particles
      for (size_t it = 0; it < par->part_maxiter; ++it) {
        for (size_t i = 0; i < sc.size(); ++i) {
          PartForce(&(xx[sx[i]]), sx[i+1] - sx[i],
                    &(ll[sl[i]]), sl[i+1] - sl[i], &(ff[sx[i]]));
        }

        // report error
        size_t dr = std::max<size_t>(1, 
            par->part_maxiter / par->part_report_fr);
        if (m.IsRoot() && (it % dr == 0 || it + 1 == par->part_maxiter)) {
          Vect h = GetCellSize();
          Scal hm = h.norminf();
          Scal tmax = 0.;
          Scal anmax = 0.;
          Scal anavg = 0; // average difference from mean angle
          Scal lmax = 0.; // maximum error in sement length
          size_t anavgn = 0;
          for (size_t i = 0; i < sc.size(); ++i) {
            // maximum force
            for (size_t j = sx[i]; j < sx[i+1]; ++j) {
              tmax = std::max(tmax, ff[j].norm());
            }
            Scal anm = 0.;
            size_t anmn = 0;
            // mean angle in string
            for (size_t j = sx[i] + 1; j + 1 < sx[i+1]; ++j) {
              anm += PS::GetAn(xx.data(), j);
              ++anmn;
            }
            anm /= anmn;
            // error in angle
            for (size_t j = sx[i] + 1; j + 1 < sx[i+1]; ++j) {
              Scal e = std::abs(anm - PS::GetAn(xx.data(), j));
              anavg += e;
              ++anavgn;
              anmax = std::max(anmax, e);
            }
            // maximum error in segment length
            for (size_t j = sx[i]; j + 1 < sx[i+1]; ++j) {
              lmax = std::max(
                  lmax, std::abs(xx[j + 1].dist(xx[j]) - par->part_h * hm));
            }
          }
          anavg /= anavgn;
          std::cout << std::setprecision(10)
              << "it=" << it 
              << " dxmax=" << tmax 
              << " anmax=" << anmax 
              << " anavg=" << anavg 
              << " lmax=" << lmax / (par->part_h * hm)
              << std::endl;
        }

        // advance
        for (size_t i = 0; i < xx.size(); ++i) {
          xx[i] += ff[i];
        }

        // copy back to field
        if (dm || it + 1 == par->part_maxiter) {
          fcps_.Reinit(m, 0);

          for (size_t i = 0; i < sc.size(); ++i) {
            IdxCell c = sc[i];
            Vect mx = smx[i];
            Vect my = smy[i];
            Vect mc = smc[i];
            for (size_t j = sx[i]; j < sx[i+1]; ++j) {
              auto q = xx[j];
              fcp_[c][fcps_[c]++] = mc + mx * q[0] + my * q[1];
            }
          }
        }

        if (dm) {
          //// XXX: disable iter dump to avoid suspender loop
          //DumpParticles(it + 1); 
        }
      }

      // compute curvature on strings
      for (size_t i = 0; i < sc.size(); ++i) {
        sk[i] = PS::PartK(&(xx[sx[i]]), sx[i + 1] - sx[i]);
      }

      // compute curvature in cells
      fckp_.Reinit(m, 0.);
      {
        size_t i = 0;
        while (i < sc.size()) {
          IdxCell c = sc[i];
          Scal k = sk[i];
          ++i;
          size_t nk = 1;
          // average over all strings in c
          while (i < sc.size() && sc[i] == c) {
            k += sk[i];
            ++i;
            ++nk;
          }
          if (par->dim == 3) {
            k *= 2.;
          }
          fckp_[c] = k / nk;
        }
      }
      m.Comm(&fckp_);

      // compute normal
      if (par->part_n) {
        for (auto c : m.Cells()) {
          if (fcps_[c]) {
            int i = 2;
            int im = i - 1;
            int ip = i + 1;
            Vect x = fcp_[c][i];
            Vect xm = fcp_[c][im];
            Vect xp = fcp_[c][ip];
            Vect dm = x - xm;
            Vect dp = xp - x;
            Scal lm = dm.norm();
            Scal lp = dp.norm();
            // normals to segments, <nm,dm> positively oriented
            Vect nm = Vect(dm[1], -dm[0], 0.) / lm; 
            Vect np = Vect(dp[1], -dp[0], 0.) / lp;
            fc_n_[c] = (nm + np) * (-0.5);
          }
        }
      }
    }

    if (sem.Nested("part-dumplast")) {
      if (dm) {
        DumpParticles(par->part_maxiter - 1);
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
      fci_.Reinit(m, false);
      for (auto c : m.AllCells()) {
        Scal u = uc[c];
        if (u > 0. && u < 1.) {
          fci_[c] = true;
        }
      }
      for (auto f : m.SuFaces()) {
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        if (uc[cm] != uc[cp]) {
          fci_[cm] = true;
          fci_[cp] = true;
        }
      }
      // Compute normal and curvature [s]
      CalcNormal(uc, fci_, fc_n_, fck_);
      auto h = GetCellSize();
      // Reconstruct interface [s]
      for (auto c : m.SuCells()) {
        fc_a_[c] = R::GetLineA(fc_n_[c], uc[c], h);
      }
    }

    // Correction with normal from particles
    if (par->part && par->part_n) {
      Part(uc, sem);
      if (sem("parta")) {
        auto h = GetCellSize();
        for (auto c : m.AllCells()) {
          fc_a_[c] = R::GetLineA(fc_n_[c], uc[c], h);
        }
      }
    }
  }
  // Print column of datafield
  void Print(const FieldFace<Scal>& ff, std::string name) {
    using MIdx = typename M::MIdx;
    auto ibc = m.GetInBlockCells();
    auto bc = m.GetBlockCells();
    auto bf = m.GetBlockFaces();
    MIdx wb = ibc.GetBegin();
    MIdx we = ibc.GetEnd();
    MIdx wc = (wb + we - MIdx(1)) / 2;
    std::cerr << std::setw(5) << name << " = ";
    for (int i = wb[1]; i <= we[1]; ++i) {
      MIdx w = wc;
      w[1] = i;
      using Dir = typename M::Dir;
      IdxFace f = bf.GetIdx(w, Dir::j);
      std::cerr << std::setw(10) << ff[f] << " ";
    }
    std::cerr << std::endl;
  }
  // Print column of datafield
  void Print(const FieldCell<Scal>& fc, std::string name) {
    using MIdx = typename M::MIdx;
    auto ibc = m.GetInBlockCells();
    auto bc = m.GetBlockCells();
    MIdx wb = ibc.GetBegin();
    MIdx we = ibc.GetEnd();
    MIdx wc = (wb + we - MIdx(1)) / 2;
    std::cerr << std::setw(5) << name << " = ";
    for (int i = wb[1]; i < we[1]; ++i) {
      MIdx w = wc;
      w[1] = i;
      IdxCell c = bc.GetIdx(w);
      std::cerr << std::setw(10) << fc[c] << " ";
    }
    std::cerr << std::endl;
  }

  void MakeIteration() override {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      auto& uc = fc_u_.iter_curr;
      const Scal dt = this->GetTimeStep();
      for (auto c : m.Cells()) {
        uc[c] = fc_u_.time_prev[c] +  // previous time step
            dt * (*fcs_)[c]; // source
      }
    }

    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    auto& bf = m.GetBlockFaces();
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
      // TODO: fluxes computed twice, consider buffer
      if (sem("adv")) {
        size_t d = dd[id]; // direction as index
        Dir md(d); // direction as Dir
        MIdx wd(md); // offset in direction d
        auto& uc = fc_u_.iter_curr;
        auto& bc = m.GetBlockCells();
        auto& bf = m.GetBlockFaces();
        auto h = GetCellSize();
        auto& ffv = *ffv_; // [f]ield [f]ace [v]olume flux
        const Scal dt = this->GetTimeStep();

        ffvu_.Reinit(m);

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
              ffvu_[f] = R::GetLineFlux(fc_n_[cu], fc_a_[cu], h, v, dt, d);
            } else { // Lagrange Explicit
              // upwind face
              IdxFace fu = bf.GetIdx(v > 0. ? wf - wd : wf + wd, md);
              // upwind mixture flux
              Scal vu = ffv[fu] * vsc;
              // phase 2 flux
              ffvu_[f] = R::GetLineFluxStr(fc_n_[cu], fc_a_[cu], h, v, vu, dt, d);
            }
          } else {
            ffvu_[f] = v * uc[cu];
          }
        }


        ffu_.Reinit(m);
        // interpolate field value to boundaries
        InterpolateB(uc, mfc_, ffu_, m);

        // override boundary upwind flux
        for (const auto& it : mfc_) {
          IdxFace f = it.GetIdx();
          CondFace* cb = it.GetValue().get(); 
          Scal v = ffv[f];
          if ((cb->GetNci() == 0) != (v > 0.)) {
            ffvu_[f] = v * ffu_[f];
            // XXX: adhoc
            // Alternating mul correction of flux
            // (done for bubble detachment)
            if (par->bcc_t0 > 0. && par->bcc_t1 > 0.) {
              Scal k0 = par->bcc_k0;
              Scal k1 = par->bcc_k1;
              Scal t0 = par->bcc_t0;
              Scal t1 = par->bcc_t1;
              Scal t = this->GetTime();
              Scal ts = t0 + t1;
              Scal ph = t / ts; 
              ph = ph - int(ph);
              ph *= ts;
              ffvu_[f] *= (ph < t0 ? k0 : k1);
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
          Scal qm = ffvu_[fm];
          Scal qp = ffvu_[fp];
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
          }
        }
        m.Comm(&uc);
      }
      if (sem.Nested("reconst")) {
        Rec(fc_u_.iter_curr);
      }
    }

    // Curvature from gradient of volume fraction
    if (par->curvgrad && sem("curv")) {
      ffu_ = Interpolate(fc_u_.iter_curr, mfc_, m); // [s]
      fcg_ = Gradient(ffu_, m); // [s]
      ffg_ = Interpolate(fcg_, mfvz_, m); // [i]

      fck_.Reinit(m, GetNan<Scal>()); // curvature [i]
      for (auto c : m.Cells()) {
        if (!fci_[c]) {
          continue;
        }
        Scal s = 0.;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          auto& g = ffg_[f];
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
      bool dm = par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                              this->GetTimeStep());
      if (dm && sem("dumppoly")) {
        DumpPoly();
      }
    }

    if (par->part) {
      Part(fc_u_.iter_curr, sem);
    }

    if (sem("stat")) {
      this->IncIter();
      ++count_;
    }
  }
  void FinishStep() override {
    fc_u_.time_curr = fc_u_.iter_curr;
    this->IncTime();
  }
  const FieldCell<Scal>& GetField(Layers l) const override {
    return fc_u_.Get(l);
  }
  const FieldCell<Scal>& GetAlpha() const {
    return fc_a_;
  }
  const FieldCell<Vect>& GetNormal() const {
    return fc_n_;
  }
  const FieldCell<Scal>& GetCurv() const override {
    return par->part_k ? fckp_ : fck_;
  }
  // curvature from height function
  const FieldCell<Scal>& GetCurvH() const {
    return fck_;
  }
  // curvature from particles
  const FieldCell<Scal>& GetCurvP() const {
    return fckp_;
  }
  using P::GetField;
};

} // namespace solver
