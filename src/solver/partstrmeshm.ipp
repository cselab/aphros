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
struct PartStrMeshM<M_>::Imp {
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;
  using R = Reconst<Scal>;
  static constexpr Scal kClNone = -1;

  Imp(M& m, std::shared_ptr<Par> par, const GRange<size_t>& layers) 
      : m(m), par(par), layers(layers), vfckp_(layers.size())
  {
    // particle strings
    partstr_ = std::unique_ptr<PS>(new PS(par->ps));
    vfckp_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
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

  void Seed(const Multi<const FieldCell<Scal>*>& vfcu,
            const Multi<const FieldCell<Scal>*>& vfca,
            const Multi<const FieldCell<Vect>*>& vfcn,
            const Multi<const FieldCell<bool>*>& vfci,
            const Multi<const FieldCell<Scal>*>& vfccl,
            const FieldCell<Scal>* fck) {
    using MIdx = typename M::MIdx;
    auto& bc = m.GetIndexCells();

    // clear string list
    partstr_->Clear();
    vsc_.clear();
    vsl_.clear();
    vsan_.clear();
    bool nocl = (vfcu.size() == 1 && !vfccl[0]); // no color provided
    Vect h = m.GetCellSize();

    // Seed strings in cells with interface.
    for (size_t l = 0; l < vfcu.size(); ++l) {
      auto& fcu = *vfcu[l];
      auto& fca = *vfca[l];
      auto& fcn = *vfcn[l];
      auto& fci = *vfci[l];
      auto& fccl = *vfccl[l];
      for (auto c : m.Cells()) {
        Vect xc = m.GetCenter(c);
        const Scal th = par->intth;
        if (fci[c] && fcu[c] >= th && fcu[c] <= 1. - th && 
            (nocl || fccl[c] != kClNone) &&
            (!fck || IsNan((*fck)[c]) ||
             std::abs((*fck)[c]) > 1. / (par->maxr * h[0]))) {

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
              for (size_t j = 0; j < vfcu.size(); ++j) {
                auto& fcu2 = *vfcu[j];
                auto& fca2 = *vfca[j];
                auto& fcn2 = *vfcn[j];
                auto& fci2 = *vfci[j];
                auto& fccl2 = *vfccl[j];
                if (fci2[cc] && (nocl || fccl2[cc] == fccl[c])) { // XXX
                  AppendInterface(v, m.GetCenter(cc), fcu2[cc],
                                  fca2[cc], fcn2[cc], fci2[cc], lx, ls);
                }
              }
            }

            // add string 
            partstr_->Add(Vect2(0.), Vect2(1., 0.), lx, ls);
            vsc_.push_back(c);
            vsl_.push_back(l);
            vsan_.push_back(an);
            assert(vsc_.size() == partstr_->GetNumStr());
            assert(vsl_.size() == partstr_->GetNumStr());
            assert(vsan_.size() == partstr_->GetNumStr());
          }
        }
      }
    }
  }
  void Part(const Multi<const FieldCell<Scal>*>& vfcu,
            const Multi<const FieldCell<Scal>*>& vfca,
            const Multi<const FieldCell<Vect>*>& vfcn,
            const Multi<const FieldCell<bool>*>& vfci,
            const Multi<const FieldCell<Scal>*>& vfccl,
            const FieldCell<Scal>* fck,
            const MapFace<std::shared_ptr<CondFace>>& mfc) {
    auto sem = m.GetSem("part");
    (void) mfc;

    if (sem("part-run")) {
      Seed(vfcu, vfca, vfcn, vfci, vfccl, fck);
      partstr_->Run(par->tol, par->itermax, m.IsRoot() ? par->verb : 0);
      // compute curvature
      vfckp_.resize(vfcu.size());
      for (auto& fckp : vfckp_.data()) {
        fckp.Reinit(m, GetNan<Scal>());
      }
      if (fck) {
        vfckp_[0] = *fck;
      }
      // XXX: assume strings from same cell contiguous
      auto ns = partstr_->GetNumStr();
      size_t s = 0;
      // traverse strings
      while (s < ns) {
        size_t nsc = 0; // number of strings in current cell
        IdxCell c = vsc_[s];
        size_t l = vsl_[s];
        auto& k = vfckp_[l][c];
        k = 0.; // reset sum
        while (s < ns && c == vsc_[s] && l == vsl_[s]) {
          k += partstr_->GetCurv(s);
          ++nsc;
          ++s;
        }
        k /= nsc;
        if (par->dim == 3) {
          k *= 2.;
        }
      }
      for (auto& fckp : vfckp_.data()) {
        m.Comm(&fckp);
      }
    }
  }

  // Dump particles to csv.
  // fca: plane constant
  // fcn: normal
  // n: frame index
  // t: time
  void DumpParticles(const Multi<const FieldCell<Scal>*>& vfca,
                     const Multi<const FieldCell<Vect>*>& vfcn,
                     size_t id, Scal t) {
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
          size_t l = vsl_[s];
          auto v = GetPlaneBasis(m.GetCenter(c),
                                 (*vfcn[l])[c], (*vfca[l])[c], vsan_[s]);

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
              << "dump" << " t=" << t << " to " << s << std::endl;
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
  void DumpPartInter(const Multi<const FieldCell<Scal>*>& vfca,
                     const Multi<const FieldCell<Vect>*>& vfcn,
                     size_t id, Scal t) {
    auto sem = m.GetSem("dumppartinter");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();

      for (size_t s = 0; s < partstr_->GetNumStr(); ++s) {
        // cell containing string
        IdxCell c = vsc_[s];
        size_t l = vsl_[s];
        auto v = GetPlaneBasis(m.GetCenter(c),
                               (*vfcn[l])[c], (*vfca[l])[c], vsan_[s]);
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
            << "dump" << " t=" << t << " to " << fn << std::endl;
        WriteVtkPoly<Vect>(fn, dl_, nullptr, {&dlc_}, {"c"}, 
            "Lines of interface around particles", false, 
            par->vtkbin, par->vtkmerge);
      }
    }
  }

  const FieldCell<Scal>& GetCurv(size_t l) {
    return vfckp_[l];
  }

 private:
  M& m;
  std::shared_ptr<Par> par;
  const GRange<size_t> layers;
  Multi<FieldCell<Scal>> vfckp_; // curvature from particles

  std::unique_ptr<PS> partstr_; // particle strings
  std::vector<IdxCell> vsc_; // vsc_[s] is cell of string s
  std::vector<size_t> vsl_; // vsl_[s] is layer of string s
  std::vector<Scal> vsan_; // vsan_[s] is angle of tangent (see GetPlaneBasis)

  std::vector<Vect> dpx_; // dump particles x
  std::vector<size_t> dpc_; // dump particles cell
  std::vector<Scal> dpk_; // dump particles curvature

  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<Scal> dlc_; // dump poly
};

template <class M_>
PartStrMeshM<M_>::PartStrMeshM(M& m, std::shared_ptr<Par> par, 
                               const GRange<size_t>& layers) 
    : imp(new Imp(m, par, layers)) {}

template <class M_>
PartStrMeshM<M_>::~PartStrMeshM() = default;

template <class M_>
void PartStrMeshM<M_>::Part(
    const Multi<const FieldCell<Scal>*>& vfcu,
    const Multi<const FieldCell<Scal>*>& vfca,
    const Multi<const FieldCell<Vect>*>& vfcn,
    const Multi<const FieldCell<bool>*>& vfci,
    const Multi<const FieldCell<Scal>*>& vfccl,
    const FieldCell<Scal>* fck,
    const MapFace<std::shared_ptr<CondFace>>& mfc) {
  imp->Part(vfcu, vfca, vfcn, vfci, vfccl, fck, mfc);
}

template <class M_>
void PartStrMeshM<M_>::DumpPartInter(
    const Multi<const FieldCell<Scal>*>& vfca,
    const Multi<const FieldCell<Vect>*>& vfcn,
    size_t id, Scal t) {
  imp->DumpPartInter(vfca, vfcn, id, t);
}

template <class M_>
void PartStrMeshM<M_>::DumpParticles(
    const Multi<const FieldCell<Scal>*>& vfca,
    const Multi<const FieldCell<Vect>*>& vfcn,
    size_t id, Scal t) {
  imp->DumpParticles(vfca, vfcn, id, t);
}

template <class M_>
auto PartStrMeshM<M_>::GetCurv(size_t l) -> const FieldCell<Scal>& {
  return imp->GetCurv(l);
}

} // namespace solver
