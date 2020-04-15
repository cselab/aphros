// Created by Petr Karnakov on 28.08.2019
// Copyright 2019 ETH Zurich

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <sstream>

#include "debug/isnan.h"
#include "dump/vtk.h"
#include "geom/mesh.h"
#include "reconst.h"
#include "solver.h"

#include "partstrmeshm.h"

template <class T>
void Reorder(std::vector<T>& v, const std::vector<size_t> idx) {
  std::vector<T> t = v;
  for (size_t i = 0; i < v.size(); ++i) {
    v[i] = t[idx[i]];
  }
}

template <class M_>
struct PartStrMeshM<M_>::Imp {
  static constexpr size_t dim = M::dim;
  using Vect2 = generic::Vect<Scal, 2>;
  using R = Reconst<Scal>;
  static constexpr Scal kClNone = -1;

  Imp(M& m, Par par, const GRange<size_t>& layers)
      : m(m), par(par), layers(layers), vfckp_(layers.size()) {
    par.ps.hc = m.GetCellSize().norminf();
    partstr_.reset(new PartStr<Scal>(par.ps));
    vfckp_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
  }

  // Clips parameterized line segment.
  // parametrization of segment:
  // x(s) = xc + t * s
  // s0,s1: endpoints
  // xmin,xmax: clipping limits, xmin < xmax
  // Returns true clipped segment is non-empty and updates s0 and s1.
  static bool ClipSegment(
      Scal xc, Scal t, Scal xmin, Scal xmax, Scal& s0, Scal& s1) {
    if (s0 > s1) {
      std::swap(s0, s1);
    }
    const Scal x0 = xc + t * s0;
    const Scal x1 = xc + t * s1;
    if (t >= 0) { // x0 < x1
      if (x1 <= xmin || x0 >= xmax) {
        return false;
      }
      if (x0 < xmin) {
        s0 = (xmin - xc) / t;
      }
      if (x1 > xmax) {
        s1 = (xmax - xc) / t;
      }
      return true;
    }
    // t < 0, x1 < x0
    if (x0 <= xmin || x1 >= xmax) {
      return false;
    }
    if (x1 < xmin) {
      s1 = (xmin - xc) / t;
    }
    if (x0 > xmax) {
      s0 = (xmax - xc) / t;
    }
    return true;
  }

  // Intersection between PLIC polygon and plane.
  // xc,n,a,h: PLIC cell center, normal, plane constant, cell size
  // xp,np: plane origin, normal
  // Output:
  // x0,x1: endpoints of intersection
  // Returns true if intersect.
  static bool GetInterPoly(
      Vect xc, Vect n, Scal a, Vect h, Vect xp, Vect np, Vect& x0, Vect& x1) {
    const Vect hh = h * 0.5;
    // plane equation: np.dot(x - xp) = 0
    // PLIC plane equation: n.dot(x - xc) = a

    // plane and cell do not intersect
    if ((np * hh).norm1() < std::abs(np.dot(xc - xp))) {
      return false;
    }
    // intersection point witn parameters u,v:
    //   x = xc + n * u + np * v
    // intersection conditions:
    //   n.dot(x - xc) = a
    //   np.dot(x - xc) = np.dot(xp - xc)
    // intersection conditions for parameters u,v:
    //   n.dot(n) * u + n.dot(np) * v = a
    //   np.dot(n) * u + np.dot(np) * v = np.dot(xp - xc)
    // matrix:
    const Scal m11 = n.dot(n);
    const Scal m12 = n.dot(np);
    const Scal m22 = np.dot(np);
    // determinant:
    const Scal det = m11 * m22 - sqr(m12);
    if (det == 0) {
      return false;
    }
    // rhs:
    const Scal r1 = a;
    const Scal r2 = np.dot(xp - xc);
    // solution:
    const Scal u = (r1 * m22 - r2 * m12) / det;
    const Scal v = (r2 * m11 - r1 * m12) / det;
    // point on intersection line:
    const Vect xl = xc + n * u + np * v;
    const int dim = 3;
    // intersection vector
    Vect t = n.cross(np);
    if (t.norm1() == 0) {
      return false;
    }
    t /= t.norm1();
    // line parametrization:
    //   x = xl + t * s
    Scal ss[2];
    const int im = t.abs().argmax();
    ss[0] = (xc[im] - hh[im] - xl[im]) / t[im];
    ss[1] = (xc[im] + hh[im] - xl[im]) / t[im];
    for (int i = 0; i < dim; ++i) {
      if (!ClipSegment(
              xl[i], t[i], xc[i] - hh[i], xc[i] + hh[i], ss[0], ss[1])) {
        return false;
      }
    }
    x0 = xl + t * ss[0];
    x1 = xl + t * ss[1];
    return true;
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
  std::array<Vect, 3> GetPlaneBasis(
      const Vect& xc, const Vect& n, Scal a, Scal an) {
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
    Vect mx = t; // unit in x, tangent
    Vect my = n / n.norm(); // unit in y, normal
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

  // Appends interface line of one cell.
  // v: output of GetPlaneBasis()
  // xc: cell center
  // a: plane constant
  // n: interface normal
  // in: interface flag
  // lx: nodes
  // ls: sizes
  // Output:
  // lx, ls: appended with interface element
  bool AppendInterface(
      const std::array<Vect, 3>& v, Vect xc, Scal a, const Vect& n,
      std::vector<Vect2>& lx, std::vector<size_t>& ls) {
    Vect h = m.GetCellSize(); // cell size
    std::array<Vect, 2> e; // ends of intersection
    Vect mc = v[0]; // plane center
    Vect mx = v[1]; // unit in x
    Vect my = v[2]; // unit in y
    Vect mn = mx.cross(my); // normal to plane

    auto xx = R::GetCutPoly(xc, n, a, h); // interface polygon
    if (R::GetInterPoly(xx, mc, mn, e)) { // intersection non-empty
      // if (GetInterPoly(xc, n, a, h, mc, mn, e[0], e[1])) {
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
    return false;
  }

  template <class EB>
  void Seed(
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn,
      const Multi<const FieldCell<bool>*>& vfci,
      const Multi<const FieldCell<Scal>*>& vfccl, const EB& eb) {
    // clear string list
    partstr_->Clear();
    vsc_.clear();
    vsl_.clear();
    vsan_.clear();
    const bool nocl = (layers.size() == 1 && !vfccl[0]); // no color provided

    // Seed strings in cells with interface.
    for (auto l : layers) {
      auto& fca = *vfca[l];
      auto& fcn = *vfcn[l];
      auto& fci = *vfci[l];
      auto& fccl = *vfccl[l];
      for (auto c : eb.Cells()) {
        if (!eb.IsRegular(c)) {
          continue;
        }
        const Vect xc = m.GetCenter(c);
        if (fci[c] && (nocl || fccl[c] != kClNone)) {
          // number of strings
          size_t ns = (par.dim == 2 ? 1 : par.ns);
          for (size_t s = 0; s < ns; ++s) {
            Scal an = s * M_PI / ns; // angle
            const auto v = GetPlaneBasis(xc, fcn[c], fca[c], an);

            // buffer for interface lines
            std::vector<Vect2> lx; // nodes
            std::vector<size_t> ls; // sizes

            auto append_contang = [&lx, &ls, &v, &fcn, &eb, this](
                                      IdxCell c, Vect2 xl, Vect nf) {
              const Scal h = m.GetCellSize()[0];
              const auto x1 = m.GetCenter(c)[1];
              //const Scal angledeg = 45. + std::abs(x1 - 1) * 90.;
              //const Scal angledeg = 45. + std::abs(x1 - 2) * 45.;
              //const Scal angledeg = 90. - (x1 - 1) * 45.;
              const Scal angledeg = 90. - (x1 - 1) * 90.;
              //const Scal angledeg = 90;
              const Scal angle = angledeg * M_PI / 180.;
              auto unit = [](const Vect& t) { return t / t.norm(); };
              const Vect tf = unit(unit(fcn[c]).orth(nf));
              const Vect ni = tf * std::sin(angle) - nf * std::cos(angle);
              Vect dir = ni.cross(v[1].cross(v[2]));
              if (dir.dot(nf) < 0) {
                dir *= -1;
              }
              const Vect xc = GetSpaceCoords(xl, v);
              auto pe0 = GetPlaneCoords(xc, v);
              auto pe1 = GetPlaneCoords(xc + dir * (10 * h), v);
              lx.push_back(pe0);
              lx.push_back(pe1);
              ls.push_back(2);
            };

            // Extract interface from neighbour cells.
            for (auto cc : m.Stencil5(c)) {
              for (auto j : layers) {
                auto& fca2 = *vfca[j];
                auto& fcn2 = *vfcn[j];
                auto& fci2 = *vfci[j];
                auto& fccl2 = *vfccl[j];
                if (fci2[cc] && (nocl || fccl2[cc] == fccl[c])) {
                  std::vector<Vect2> llx; // nodes
                  std::vector<size_t> lls; // sizes
                  if (AppendInterface(
                        v, m.GetCenter(cc), fca2[cc], fcn2[cc], llx, lls)) {
                    bool nearcut = false;
                    Vect nf;
                    if (c == cc) {
                      for (auto cn : eb.Stencil(cc)) {
                        if (eb.IsCut(cn)) {
                          nearcut = true;
                          nf = eb.GetNormal(cn);
                          break;
                        }
                      }
                    }
                    if (nearcut && c == cc) {
                      append_contang(cc, (llx[0] + llx[1]) * 0.5, nf);
                    } else if (!eb.IsRegular(cc)) {
                      // nop
                    } else {
                      lx.push_back(llx[0]);
                      lx.push_back(llx[1]);
                      ls.push_back(lls[0]);
                    }
                  }
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
  template <class EB>
  void Part(
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn,
      const Multi<const FieldCell<bool>*>& vfci,
      const Multi<const FieldCell<Scal>*>& vfccl, const EB& eb) {
    auto sem = m.GetSem("part");

    if (sem("part-run")) {
      Seed(vfca, vfcn, vfci, vfccl, eb);
      partstr_->Run(par.tol, par.itermax, m.IsRoot() ? par.verb : 0);
      // compute curvature
      vfckp_.resize(layers);
      for (auto& fckp : vfckp_.data()) {
        fckp.Reinit(m, GetNan<Scal>());
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
        if (par.dim == 3) {
          k *= 2.;
        }
      }
      for (auto& fckp : vfckp_.data()) {
        m.Comm(&fckp);
      }
    }
    const bool nocl = (layers.size() == 1 && !vfccl[0]); // no color provided
    if (sem("part")) {
      // finds a valid curvature value in stencil around `c`
      // selecting cells with color `cl`
      auto findcurv = [&](IdxCell c, size_t l) {
        for (auto cn : eb.Stencil(c)) {
          for (auto ln : layers) {
            if (!IsNan(vfckp_[ln][cn]) &&
                (nocl || (*vfccl[ln])[cn] == (*vfccl[l])[c])) {
              return vfckp_[ln][cn];
            }
          }
        }
        return GetNan<Scal>();
      };
      for (auto c : eb.Cells()) {
        if (eb.IsCut(c)) {
          for (auto l : layers) {
            vfckp_[l][c] = findcurv(c, l);
          }
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
  void DumpParticles(
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn, size_t id, Scal t) {
    auto sem = m.GetSem("partdump");
    struct {
      std::vector<Vect> dpx; // dump particles x
      std::vector<size_t> dpc; // dump particles cell
      std::vector<Scal> dpk; // dump particles curvature
    } * ctx(sem);
    auto& dpx = ctx->dpx;
    auto& dpc = ctx->dpc;
    auto& dpk = ctx->dpk;

    size_t it = 1;
    if (1) { // TODO: revise frames
      if (sem("local")) {
        dpx.clear();
        dpc.clear();
        dpk.clear();

        // loop over strings
        for (size_t s = 0; s < partstr_->GetNumStr(); ++s) {
          // cell
          IdxCell c = vsc_[s];
          size_t l = vsl_[s];
          auto v = GetPlaneBasis(
              m.GetCenter(c), (*vfcn[l])[c], (*vfca[l])[c], vsan_[s]);

          auto p = partstr_->GetStr(s);
          const Vect2* xx = p.first;
          size_t sx = p.second;

          size_t ic = m.GetHash(c);
          for (size_t i = 0; i < sx; ++i) {
            auto x = GetSpaceCoords(xx[i], v);
            dpx.push_back(x);
            dpc.push_back(ic);
            dpk.push_back(partstr_->GetCurv(s));
          }
        }

        // comm
        using TV = typename M::template OpCatT<Vect>;
        using TI = typename M::template OpCatT<size_t>;
        using TS = typename M::template OpCatT<Scal>;
        m.Reduce(std::make_shared<TV>(&dpx));
        m.Reduce(std::make_shared<TI>(&dpc));
        m.Reduce(std::make_shared<TS>(&dpk));
      }
      if (sem("write")) {
        if (m.IsRoot()) {
          std::string s =
              GetDumpName("partit", ".csv", id, par.dump_fr > 1 ? it : -1);
          std::cout << std::fixed << std::setprecision(8) << "dump"
                    << " t=" << t << " to " << s << std::endl;
          std::ofstream o;
          o.open(s);
          o.precision(20);
          o << "x,y,z,c,k\n";

          std::vector<size_t> idx(dpc.size());
          std::iota(idx.begin(), idx.end(), 0);
          std::stable_sort(
              idx.begin(), idx.end(),
              [&dpc](size_t i1, size_t i2) { return dpc[i1] < dpc[i2]; });

          for (size_t j = 0; j < dpx.size(); ++j) {
            const size_t i = idx[j];
            Vect x = dpx[i];
            o << x[0] << "," << x[1] << "," << x[2] << "," << dpc[i] << ","
              << dpk[i] << "\n";
          }
        }
      }
    }
  }
  // Dumps interface around particle strings
  void DumpPartInter(
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn, size_t id, Scal t) {
    auto sem = m.GetSem("dumppartinter");
    struct {
      std::vector<std::vector<Vect>> dl; // lines
      std::vector<Scal> dlc; // cell indices
    } * ctx(sem);
    auto& dl = ctx->dl;
    auto& dlc = ctx->dlc;

    if (sem("local")) {
      for (size_t s = 0; s < partstr_->GetNumStr(); ++s) {
        // cell containing string
        IdxCell c = vsc_[s];
        size_t l = vsl_[s];
        auto v = GetPlaneBasis(
            m.GetCenter(c), (*vfcn[l])[c], (*vfca[l])[c], vsan_[s]);
        auto in = partstr_->GetInter(s);
        size_t i = 0;
        size_t hc = m.GetHash(c);
        for (size_t l = 0; l < in.n; ++l) {
          std::vector<Vect> xx;
          for (size_t k = 0; k < in.s[l]; ++k) {
            xx.push_back(GetSpaceCoords(in.x[i], v));
            ++i;
          }
          dl.push_back(xx);
          dlc.push_back(hc);
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::vector<size_t> idx(dlc.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {
          return dlc[i1] < dlc[i2];
        });

        Reorder(dl, idx);
        Reorder(dlc, idx);

        std::string fn = GetDumpName("sp", ".vtk", id);
        std::cout << std::fixed << std::setprecision(8) << "dump"
                  << " t=" << t << " to " << fn << std::endl;
        WriteVtkPoly<Vect>(
            fn, dl, nullptr, {&dlc}, {"c"},
            "Lines of interface around particles", false, par.vtkbin,
            par.vtkmerge);
      }
    }
  }

  Multi<const FieldCell<Scal>*> GetCurv() {
    return vfckp_;
  }

 private:
  M& m;
  Par par;
  const GRange<size_t> layers;
  Multi<FieldCell<Scal>> vfckp_; // curvature from particles

  std::unique_ptr<PartStr<Scal>> partstr_; // particle strings
  std::vector<IdxCell> vsc_; // vsc_[s] is cell of string s
  std::vector<size_t> vsl_; // vsl_[s] is layer of string s
  std::vector<Scal> vsan_; // vsan_[s] is angle of tangent (see GetPlaneBasis)
};

template <class M_>
PartStrMeshM<M_>::PartStrMeshM(M& m, Par par, const GRange<size_t>& layers)
    : imp(new Imp(m, par, layers)) {}

template <class M_>
PartStrMeshM<M_>::~PartStrMeshM() = default;

template <class M_>
template <class EB>
void PartStrMeshM<M_>::Part(
    const Multi<const FieldCell<Scal>*>& vfca,
    const Multi<const FieldCell<Vect>*>& vfcn,
    const Multi<const FieldCell<bool>*>& vfci,
    const Multi<const FieldCell<Scal>*>& vfccl, const EB& eb) {
  imp->Part(vfca, vfcn, vfci, vfccl, eb);
}

template <class M_>
void PartStrMeshM<M_>::DumpPartInter(
    const Multi<const FieldCell<Scal>*>& vfca,
    const Multi<const FieldCell<Vect>*>& vfcn, size_t id, Scal t) {
  imp->DumpPartInter(vfca, vfcn, id, t);
}

template <class M_>
void PartStrMeshM<M_>::DumpParticles(
    const Multi<const FieldCell<Scal>*>& vfca,
    const Multi<const FieldCell<Vect>*>& vfcn, size_t id, Scal t) {
  imp->DumpParticles(vfca, vfcn, id, t);
}

template <class M_>
auto PartStrMeshM<M_>::GetCurv() -> Multi<const FieldCell<Scal>*> {
  return imp->GetCurv();
}
