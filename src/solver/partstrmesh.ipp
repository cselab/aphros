#include <memory>
#include <cmath>
#include <sstream>

#include "geom/mesh.h"
#include "solver.h"
#include "reconst.h"
#include "debug/isnan.h"
#include "dump/vtk.h"

namespace solver {

template <class M_>
struct PartStrMesh<M_>::Imp {
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;
  using R = Reconst<Scal>;

  Imp(M& m, std::shared_ptr<Par> par) : m(m), par(par), fckp_(m, 0) {
    // particle strings
    partstr_ = std::unique_ptr<PS>(new PS(par->ps));
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
  // mn: unit vector to form orthonormal positively oriented basis <mx,my,mn> 
  std::array<Vect, 3> GetPlaneBasis(const Vect& xc, const Vect& n, 
                                    Scal a, Scal an) {
    Vect h = m.GetCellSize();
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
  Vect2 GetPlaneCoords(const Vect& x, const std::array<Vect, 3>& v) {
    auto mc = v[0];
    auto mx = v[1];
    auto my = v[2];
    return Vect2((x - mc).dot(mx), (x - mc).dot(my));
  }
  // Convert from plane to space coordinates.
  // (xl,yl,0): plane coordinates
  // v: output of GetPlaneBasis()
  // Returns:
  // x: space coordinates
  Vect GetSpaceCoords(const Vect2& xl, const std::array<Vect, 3>& v) {
    auto mc = v[0];
    auto mx = v[1];
    auto my = v[2];
    return mc + mx * xl[0] + my * xl[1];
  }

  // Appends interface line
  bool AppendInterfaceLine(const std::array<Vect, 3>& v,
                           Vect xc, Scal u, Scal a, const Vect& n, bool in,
                           std::vector<Vect2>& lx, std::vector<size_t>& ls) { 
    Vect h = m.GetCellSize(); // cell size
    const Scal th = par->intth;

    if (in && u >= th && u <= 1. - th) {
      auto xx = R::GetCutPoly(xc, n, a, h); // interface polygon
      std::array<Vect, 2> e; // ends of intersection 
      Vect mc = v[0];  // plane center
      Vect mx = v[1];  // unit in x
      Vect my = v[2];  // unit in y
      Vect mn = mx.cross(my); // normal to plane
      if (R::GetInterPoly(xx, mc, mn, e)) { // intersection non-empty
        // interface normal 
        auto pn = GetPlaneCoords(mc + n, v);
        // line ends
        auto pe0 = GetPlaneCoords(e[0], v);
        auto pe1 = GetPlaneCoords(e[1], v);
        // make <pn,pe1-pe0> positively oriented
        if (pn.cross_third(pe1 - pe0) < 0.) {
          std::swap(pe0, pe1);
        }
        // polygon of fluid volume
        lx.push_back(pe0);
        lx.push_back(pe1);
        ls.push_back(2);
        return true;
      }
    }
    return false;
  }
  // Appends interface volume
  bool AppendInterfaceVolume(const std::array<Vect, 3>& v,
                             Vect xc, Scal u, Scal a, const Vect& n, bool in,
                             std::vector<Vect2>& lx, std::vector<size_t>& ls) { 
    Vect h = m.GetCellSize(); // cell size

    // XXX: adhoc 2d
    // cell contour polygon // TODO: general for 3d cell
    std::vector<Vect2> pc;
    pc.push_back(GetPlaneCoords(xc + Vect(-h[0], -h[1], 0.) * 0.5, v));
    pc.push_back(GetPlaneCoords(xc + Vect(h[0], -h[1], 0.) * 0.5, v));
    pc.push_back(GetPlaneCoords(xc + Vect(h[0], h[1], 0.) * 0.5, v));
    pc.push_back(GetPlaneCoords(xc + Vect(-h[0], h[1], 0.) * 0.5, v));

    const Scal th = par->intth;

    if (in && u >= th && u <= 1. - th) {
      auto xx = R::GetCutPoly(xc, n, a, h); // interface polygon
      std::array<Vect, 2> e; // ends of intersection 
      Vect mc = v[0];  // plane center
      Vect mx = v[1];  // unit in x
      Vect my = v[2];  // unit in y
      Vect mn = mx.cross(my); // normal to plane
      if (R::GetInterPoly(xx, mc, mn, e)) { // intersection non-empty
        // interface normal 
        auto pn = GetPlaneCoords(mc + n, v);
        // line ends
        auto pe0 = GetPlaneCoords(e[0], v);
        auto pe1 = GetPlaneCoords(e[1], v);
        // make <pn,pe1-pe0> positively oriented
        if (pn.cross_third(pe1 - pe0) < 0.) {
          std::swap(pe0, pe1);
        }
        // polygon of fluid volume
        auto pv = R::GetCutPoly(pc, {pe0, pe1});
        lx.insert(lx.end(), pv.begin(), pv.end());
        ls.push_back(pv.size());
        return true;
      } else {
        lx.insert(lx.end(), pc.begin(), pc.end());
        ls.push_back(pc.size());
        return true;
      }
    } else if (u > 1. - th) {
      lx.insert(lx.end(), pc.begin(), pc.end());
      ls.push_back(pc.size());
        return true;
    } 
    return false;
  }
  // Appends interface element (line or volume) of one cell.
  // v: output of GetPlaneBasis()
  // xc: cell center
  // u: volume fraction
  // a: plane constant
  // n: interface normal
  // in: interface flag
  // lx: nodes
  // ls: sizes
  // Output:
  // lx, ls: appended with interface element
  bool AppendInterface(const std::array<Vect, 3>& v,
                       Vect xc, Scal u, Scal a, const Vect& n, bool in,
                       std::vector<Vect2>& lx, std::vector<size_t>& ls) { 
    switch (par->attrreconst) { 
      case AR::line:
        return AppendInterfaceLine(v, xc, u, a, n, in, lx, ls);
      case AR::volume:
        return AppendInterfaceVolume(v, xc, u, a, n, in, lx, ls);
      default:
        throw std::runtime_error("AppendInterface(): Unknown attrreconst");
    }
  }

  void Seed(const FieldCell<Scal>& fcu, const FieldCell<Scal>& fca, 
             const FieldCell<Vect>& fcn, const FieldCell<bool>& fci) {
    using MIdx = typename M::MIdx;
    auto& bc = m.GetIndexCells();

    // clear string list
    partstr_->Clear();
    vsc_.clear();
    vsan_.clear();

    // Seed strings in cells with interface.
    for (auto c : m.Cells()) {
      Vect xc = m.GetCenter(c);
      const Scal th = par->intth;
      if (fci[c] && fcu[c] >= th && fcu[c] <= 1. - th) {
        // number of strings
        size_t ns = (par->dim == 2 ? 1 : par->ns);
        for (size_t s = 0; s < ns; ++s) {
          Scal an = s * M_PI / ns; // angle
          auto v = GetPlaneBasis(xc, fcn[c], fca[c], an);

          // block of offsets to neighbours in stencil [-sw,sw]
          const int sw = 2; // stencil halfwidth
          const int sn = sw * 2 + 1; // stencil size
          GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, par->dim == 2 ? 0 : -sw), 
                                  MIdx(sn, sn, par->dim == 2 ? 1 : sn)); 

          auto w = bc.GetMIdx(c);
          // buffer for interface lines
          std::vector<Vect2> lx; // nodes
          std::vector<size_t> ls; // sizes
          // Extract interface from neighbour cells.
          for (auto wo : bo) {
            IdxCell cc = bc.GetIdx(w + wo); // neighbour cell
            AppendInterface(v, m.GetCenter(cc), fcu[cc], fca[cc], fcn[cc],
                            fci[cc], lx, ls);
          }

          // add string 
          partstr_->Add(Vect2(0.), Vect2(1., 0.), lx, ls);
          vsc_.push_back(c);
          vsan_.push_back(an);
          assert(vsc_.size() == partstr_->GetNumStr());
          assert(vsan_.size() == partstr_->GetNumStr());
        }
      }
    }
  }
  void Part(const FieldCell<Scal>& uc, 
      FieldCell<Scal>& fca, FieldCell<Vect>& fcn, FieldCell<bool>& fci,
      const MapFace<std::shared_ptr<CondFace>>& mfc) {
    auto sem = m.GetSem("part");

    if (sem("part-comm")) {
      m.Comm(&fca);
      m.Comm(&fcn);
    }

    if (sem("part-run")) {
      // XXX: adhoc
      // reflection at boundaries
      if (par->bcc_reflect) {
        BcReflect(fca, mfc, m);
        BcReflect(fcn, mfc, m);
      }
      Seed(uc, fca, fcn, fci);
      partstr_->Run(par->tol, par->itermax,
                    m.IsRoot() ? par->verb : 0);
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
        if (fci[c] && IsNan(fckp_[c])) {
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
      if (nan && par->verb) {
        std::stringstream s;
        s.precision(16);
        s << "nan curvature in all neighbours of " << nan 
            << " cells, one at x=" << m.GetCenter(cnan)
            << " with vf=" << uc[cnan];
        std::cout << s.str() << std::endl;
      }
      m.Comm(&fckp_);
    }
  }

  // Dump particles to csv.
  // fca: plane constant
  // fcn: normal
  // n: frame index
  // t: time
  // dt: timestep
  void DumpParticles(FieldCell<Scal>& fca, FieldCell<Vect>& fcn,
                     size_t id, Scal t, Scal dt) {
    auto sem = m.GetSem("partdump");
    size_t it = 1;
    if (1) { // TODO: revise frames
      if (sem("local")) {
        dpx_.clear();
        dpc_.clear();
        dpk_.clear();

        // loop over strings
        for (size_t s = 0; s < partstr_->GetNumStr(); ++s) {
          // cell
          IdxCell c = vsc_[s];
          auto v = GetPlaneBasis(m.GetCenter(c), fcn[c], fca[c], vsan_[s]);

          auto p = partstr_->GetStr(s);
          const Vect2* xx = p.first;
          size_t sx = p.second;

          size_t ic = m.GetHash(c);
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
          std::string s = GetDumpName("partit", ".csv", id,
                                      par->dump_fr > 1 ? it : -1);
          std::cout << std::fixed << std::setprecision(8)
              << "dump" 
              << " t=" << t + dt
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
  // Dumps interface around particle strings
  void DumpPartInter(FieldCell<Scal>& fca, FieldCell<Vect>& fcn,
                     size_t id, Scal t, Scal dt) {
    auto sem = m.GetSem("dumppartinter");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();

      for (size_t s = 0; s < partstr_->GetNumStr(); ++s) {
        // cell containing string
        IdxCell c = vsc_[s];
        auto v = GetPlaneBasis(m.GetCenter(c), fcn[c], fca[c], vsan_[s]);
        auto in = partstr_->GetInter(s);
        size_t i = 0;
        size_t hc = m.GetHash(c);
        for (size_t l = 0; l < in.n; ++l) {
          std::vector<Vect> xx;
          for (size_t k = 0; k < in.s[l]; ++k) {
            xx.push_back(GetSpaceCoords(in.x[i], v));
            ++i;
          }
          dl_.push_back(xx);
          dlc_.push_back(hc);
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string fn = GetDumpName("sp", ".vtk", id);
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " t=" << t + dt
            << " to " << fn << std::endl;
        WriteVtkPoly(fn, dl_, {&dlc_}, {"c"}, 
            "Lines of interface around particles");
      }
    }
  }

  const FieldCell<Scal>& GetCurv() {
    return fckp_;
  }

 private:
  M& m;
  std::shared_ptr<Par> par;
  FieldCell<Scal> fckp_; // curvature from particles

  std::unique_ptr<PS> partstr_; // particle strings
  std::vector<IdxCell> vsc_; // vsc_[s] is cell of string s
  std::vector<Scal> vsan_; // vsan_[s] is angle of tangent (see GetPlaneBasis)

  std::vector<Vect> dpx_; // dump particles x
  std::vector<size_t> dpc_; // dump particles cell
  std::vector<Scal> dpk_; // dump particles curvature

  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<Scal> dlc_; // dump poly
};

template <class M_>
PartStrMesh<M_>::PartStrMesh(M& m, std::shared_ptr<Par> par) 
    : imp(new Imp(m, par)) {}

template <class M_>
PartStrMesh<M_>::~PartStrMesh() = default;

template <class M_>
void PartStrMesh<M_>::Part(const FieldCell<Scal>& uc, 
    FieldCell<Scal>& fca, FieldCell<Vect>& fcn, FieldCell<bool>& fci,
    const MapFace<std::shared_ptr<CondFace>>& mfc) {
  imp->Part(uc, fca, fcn, fci, mfc);
}

template <class M_>
void PartStrMesh<M_>::DumpPartInter(
    FieldCell<Scal>& fca, FieldCell<Vect>& fcn, size_t id, Scal t, Scal dt) {
  imp->DumpPartInter(fca, fcn, id, t, dt);
}

template <class M_>
void PartStrMesh<M_>::DumpParticles(
    FieldCell<Scal>& fca, FieldCell<Vect>& fcn, size_t id, Scal t, Scal dt) {
  imp->DumpParticles(fca, fcn, id, t, dt);
}

template <class M_>
auto PartStrMesh<M_>::GetCurv() -> const FieldCell<Scal>& {
  return imp->GetCurv();
}

} // namespace solver
