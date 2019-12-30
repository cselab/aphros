#pragma once

#include <limits>
#include <stdexcept>

#include "dump/vtk.h"
#include "reconst.h"
#include "solver.h"

template <class T>
class FieldEmbed {
 public:
  using Value = T;

  FieldEmbed() = default;
  FieldEmbed(const FieldEmbed&) = default;
  FieldEmbed(FieldEmbed&&) = default;
  FieldEmbed(const FieldCell<Value>& fc, const FieldFace<Value>& ff)
      : dc_(fc), df_(ff) {}
  FieldEmbed& operator=(const FieldEmbed& o) = default;
  FieldEmbed& operator=(FieldEmbed&& o) = default;
  template <class M>
  explicit FieldEmbed(const M& m) : dc_(m), df_(m) {}
  template <class M>
  explicit FieldEmbed(const M& m, const Value& v) : dc_(m, v), df_(m, v) {}
  template <class M>
  void Reinit(const M& m) {
    dc_.Reinit(m);
    df_.Reinit(m);
  }
  template <class M>
  void Reinit(const M& m, const Value& v) {
    dc_.Reinit(m, v);
    df_.Reinit(m, v);
  }
  Value& operator[](IdxCell c) {
    return dc_[c];
  }
  const Value& operator[](IdxCell c) const {
    return dc_[c];
  }
  Value& operator[](IdxFace f) {
    return df_[f];
  }
  const Value& operator[](IdxFace f) const {
    return df_[f];
  }
  FieldCell<Value>& GetFieldCell() {
    return dc_;
  }
  const FieldCell<Value>& GetFieldCell() const {
    return dc_;
  }
  FieldFace<Value>& GetFieldFace() {
    return df_;
  }
  const FieldFace<Value>& GetFieldFace() const {
    return df_;
  }

 private:
  FieldCell<Value> dc_;
  FieldFace<Value> df_;
};

// Embedded boundaries.
template <class M_>
class Embed {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using R = Reconst<Scal>;
  static constexpr size_t dim = M::dim;

 public:
  enum class Type { regular, cut, excluded };
  // fnl: level-set function on nodes, interface at fnl=0
  Embed(M& m, const FieldNode<Scal>& fnl) : m(m), fnl_(fnl) {
    InitFaces(fnl_, fft_, ffpoly_, ffs_, m);
    InitCells(fnl_, ffs_, fct_, fcn_, fca_, fcs_, fcv_, m);
  }
  Type GetType(IdxCell c) const {
    return fct_[c];
  }
  Type GetType(IdxFace f) const {
    return fft_[f];
  }
  // Returns outer normal (towards excluded domain) in cut cells.
  Vect GetNormal(IdxCell c) const {
    return fcn_[c];
  }
  Vect GetNormal(IdxFace f) const {
    return m.GetNormal(f);
  }
  Scal GetAlpha(IdxCell c) const {
    return fca_[c];
  }
  Scal GetArea(IdxFace f) const {
    return ffs_[f];
  }
  Scal GetArea(IdxCell c) const {
    return fcs_[c];
  }
  Scal GetFaceOffset(IdxCell c, size_t nci) const {
    IdxFace f = m.GetFace(c, nci);
    return (GetFaceCenter(f) - GetCellCenter(c)).dot(m.GetNormal(f));
  }
  Scal GetFaceOffset(IdxCell c) const {
    return (GetFaceCenter(c) - GetCellCenter(c)).dot(GetNormal(c));
  }
  Scal GetVolume(IdxCell c) const {
    return fcv_[c];
  }
  std::vector<Vect> GetFacePoly(IdxFace f) const {
    switch (fft_[f]) {
      case Type::regular:
        return GetRegularFacePoly(f, m);
      case Type::cut:
        return ffpoly_[f];
      default:
        return std::vector<Vect>();
    }
  }
  Vect GetFaceCenter(IdxCell c) const {
    if (fct_[c] == Type::cut) {
      return GetBaryCenter(
          R::GetCutPoly(m.GetCenter(c), fcn_[c], fca_[c], m.GetCellSize()));
    }
    return GetNan<Vect>();
  }
  Vect GetFaceCenter(IdxFace f) const {
    switch (fft_[f]) {
      case Type::regular:
        return m.GetCenter(f);
      case Type::cut:
        return GetBaryCenter(ffpoly_[f]);
      default:
        return GetNan<Vect>();
    }
  }
  Vect GetCellCenter(IdxCell c) const {
    switch (fct_[c]) {
      case Type::regular:
        return m.GetCenter(c);
      case Type::cut: {
        const Scal w = GetArea(c);
        Vect sum = GetFaceCenter(c) * w;
        Scal sumw = w;
        for (auto q : m.Nci(c)) {
          auto f = m.GetFace(c, q);
          if (fft_[f] != Type::excluded) {
            const Scal w = GetArea(f);
            sum += GetFaceCenter(f) * w;
            sumw += w;
          }
        }
        return sum / sumw;
      }
      default:
        return GetNan<Vect>();
    }
  }

  // feu: field on embedded boundaries [a]
  // Returns:
  // field on cells [a]
  FieldCell<Scal> Interpolate(const FieldEmbed<Scal>& feu) const {
    FieldCell<Scal> fcu(m, 0); // FIXME should be nan
    for (auto c : m.AllCells()) {
      switch (fct_[c]) {
        case Type::regular: {
          Scal sum = 0;
          Scal sumw = 0;
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            sum += feu[f];
            sumw += 1.;
          }
          fcu[c] = sum / sumw;
          break;
        }
        case Type::cut: {
          const Scal w = 1 / std::abs(GetFaceOffset(c));
          Scal sum = feu[c] * w;
          Scal sumw = w;
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            if (fft_[f] == Type::regular || fft_[f] == Type::cut) {
              const Scal w = 1 / std::abs(GetFaceOffset(c, q));
              sum += feu[f] * w;
              sumw += w;
            }
          }
          fcu[c] = sum / sumw;
          break;
        }
        case Type::excluded:
          fcu[c] = 0;
          break;
      }
    }
    return fcu;
  }
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  FieldEmbed<Scal> Interpolate(
      const FieldCell<Scal>& fcu, size_t bc, Scal bcv) const {
    FieldEmbed<Scal> feu(m, 0); // FIXME should be nan
    for (auto f : m.SuFaces()) {
      switch (fft_[f]) {
        case Type::regular:
        case Type::cut: {
          IdxCell cm = m.GetCell(f, 0);
          IdxCell cp = m.GetCell(f, 1);
          Scal a = 0.5;
          feu[f] = fcu[cm] * (1 - a) + fcu[cp] * a;
          break;
        }
        case Type::excluded: {
          feu[f] = 0;
          break;
        }
      }
    }
    for (auto c : m.SuCells()) {
      switch (fct_[c]) {
        case Type::cut: {
          if (bc == 0) {
            feu[c] = bcv;
          } else if (bc == 1) {
            feu[c] = fcu[c] + GetFaceOffset(c) * bcv;
          } else {
            throw std::runtime_error("unknown bc=" + std::to_string(bc));
          }
          break;
        }
        case Type::regular:
        case Type::excluded: {
          feu[c] = 0;
          break;
        }
      }
    }
    return feu;
  }
  void DumpPoly() const {
    const std::string fn = GetDumpName("eb", ".vtk", 0);
    DumpPoly(fn, ffs_, fft_, fcs_, fct_, fcn_, fca_, ffpoly_, m);
  }

 private:
  // Dump cut polygons
  // ffs: face area for which f > 0
  // fft: type of faces
  // fcs: polygon area
  // fct: cell types
  // fcn: normals
  // fca: plane constants
  // ffpoly: polygon representing f < 0
  static void DumpPoly(
      const std::string fn, const FieldFace<Scal>& ffs,
      const FieldFace<Type>& fft, const FieldCell<Scal>& fcs,
      const FieldCell<Type>& fct, const FieldCell<Vect>& fcn,
      const FieldCell<Scal>& fca, const FieldFace<std::vector<Vect>>& ffpoly,
      M& m) {
    auto sem = m.GetSem("dumppoly");
    struct {
      std::vector<std::vector<Vect>> dl; // polygon
      std::vector<Scal> dld; // direction
      std::vector<Scal> dls; // area
      std::vector<Scal> dlf; // 0: cell, 1: face
    } * ctx(sem);
    auto& dl = ctx->dl;
    auto& dld = ctx->dld;
    auto& dls = ctx->dls;
    auto& dlf = ctx->dlf;
    if (sem("local")) {
      for (auto f : m.Faces()) {
        if (fct[m.GetCell(f, 0)] == Type::cut ||
            fct[m.GetCell(f, 1)] == Type::cut) {
          size_t d(m.GetIndexFaces().GetDir(f));
          if (fft[f] == Type::cut || fft[f] == Type::regular) {
            if (fft[f] == Type::cut) {
              dl.push_back(ffpoly[f]);
            } else {
              dl.push_back(GetRegularFacePoly(f, m));
            }
            dld.push_back(d);
            dls.push_back(ffs[f]);
            dlf.push_back(1);
          }
        }
      }

      for (auto c : m.Cells()) {
        if (fct[c] == Type::cut) {
          auto xx =
              R::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], m.GetCellSize());
          dl.push_back(xx);
          dld.push_back(3);
          dls.push_back(fcs[c]);
          dlf.push_back(0);
        }
      }

      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dld));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::cout << std::fixed << std::setprecision(8) << "dump"
                  << " to " << fn << std::endl;
        WriteVtkPoly<Vect>(
            fn, dl, nullptr, {&dld, &dls, &dlf}, {"dir", "area", "face"},
            "Embedded boundary", true, true, true);
      }
    }
  }
  static Scal GetTriangleArea(Vect xa, Vect xb, Vect xc) {
    return (xc - xa).cross(xb - xa).norm() * 0.5;
  }
  static Vect GetBaryCenter(const std::vector<Vect>& xx) {
    Vect sum(0);
    Scal suma = 0;
    for (size_t i = 2; i < xx.size(); ++i) {
      const Scal a = GetTriangleArea(xx[0], xx[i - 1], xx[i]);
      const Vect x = (xx[0] + xx[i - 1] + xx[i]) / 3.;
      sum += x * a;
      suma += a;
    }
    return sum / suma;
  }
  // Returns point at which interpolant has value 0.
  // x0,x1: points
  // f0,f1: values
  static Vect GetIso(Vect x0, Vect x1, Scal f0, Scal f1) {
    return (x0 * f1 - x1 * f0) / (f1 - f0);
  }
  static std::vector<Vect> GetRegularFacePoly(IdxFace f, const M& m) {
    std::vector<Vect> xx;
    for (size_t e = 0; e < m.GetNumNodes(f); ++e) {
      auto n = m.GetNode(f, e);
      xx.push_back(m.GetNode(n));
    }
    return xx;
  }
  // Determines the face types, constructs polygons and computes fractions.
  // fnl: level-set on nodes, fnl > 0 inside regular cells [a]
  // Output:
  // fft: type of faces
  // ffpoly: if fft=1, polygon representing fnl > 0; otherwise empty
  // ffs: face area for which fnl > 0
  static void InitFaces(
      const FieldNode<Scal>& fnl, FieldFace<Type>& fft,
      FieldFace<std::vector<Vect>>& ffpoly, FieldFace<Scal>& ffs, const M& m) {
    fft.Reinit(m);
    ffpoly.Reinit(m);
    ffs.Reinit(m, 0);
    for (auto f : m.AllFaces()) {
      const size_t em = m.GetNumNodes(f);
      std::vector<Vect> xx;
      bool cut = false;
      for (size_t e = 0; e < em; ++e) {
        size_t ep = (e + 1) % em;
        IdxNode n = m.GetNode(f, e);
        IdxNode np = m.GetNode(f, ep);
        Scal l = fnl[n];
        Scal lp = fnl[np];
        Vect x = m.GetNode(n);
        Vect xp = m.GetNode(np);
        if (l > 0) {
          xx.push_back(x);
        }
        if ((l < 0) != (lp < 0)) {
          xx.push_back(GetIso(x, xp, l, lp));
          cut = true;
        }
      }
      fft[f] = (cut ? Type::cut : xx.empty() ? Type::excluded : Type::regular);
      switch (fft[f]) {
        case Type::regular:
          ffs[f] = m.GetArea(f);
          break;
        case Type::cut:
          ffs[f] = std::abs(R::GetArea(xx, m.GetNormal(f)));
          ffpoly[f] = xx;
          break;
        case Type::excluded:
          ffs[f] = 0;
          break;
      }
    }
  }
  // Determines the cell types, normals and plane constants.
  // fnl: level-set on nodes [a]
  // ffs: face area for which f > 0
  // Output:
  // fct: cell types
  // fcn: normals
  // fca: plane constants
  // fcs: polygon area
  // fcv: cell volume
  static void InitCells(
      const FieldNode<Scal>& fnl, const FieldFace<Scal>& ffs,
      FieldCell<Type>& fct, FieldCell<Vect>& fcn, FieldCell<Scal>& fca,
      FieldCell<Scal>& fcs, FieldCell<Scal>& fcv, const M& m) {
    fct.Reinit(m), GetNan<Scal>();
    fcn.Reinit(m), GetNan<Scal>();
    fca.Reinit(m, GetNan<Scal>());
    fcs.Reinit(m, GetNan<Scal>());
    fcv.Reinit(m, GetNan<Scal>());
    for (auto c : m.AllCells()) {
      size_t q = 0; // number of nodes with fnl > 0
      const size_t mi = m.GetNumNodes(c);
      for (size_t i = 0; i < mi; ++i) {
        IdxNode n = m.GetNode(c, i);
        if (fnl[n] >= 0) {
          ++q;
        }
      }
      fct[c] = (q == mi ? Type::regular : q > 0 ? Type::cut : Type::excluded);

      if (fct[c] == Type::cut) {
        // calc normal
        {
          Vect n(0);
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            n += m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q);
          }
          fcn[c] = -n / n.norm();
        }

        // calc plane constant
        {
          Scal a = 0;
          Scal aw = 0;
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            const size_t em = m.GetNumNodes(f);
            for (size_t e = 0; e < em; ++e) {
              size_t ep = (e + 1) % em;
              IdxNode n = m.GetNode(f, e);
              IdxNode np = m.GetNode(f, ep);
              Scal l = fnl[n];
              Scal lp = fnl[np];
              Vect x = m.GetNode(n);
              Vect xp = m.GetNode(np);

              if ((l < 0) != (lp < 0)) {
                a += fcn[c].dot(GetIso(x, xp, l, lp) - m.GetCenter(c));
                aw += 1.;
              }
            }
          }
          fca[c] = a / aw;
        }

        const auto h = m.GetCellSize();

        auto xx = R::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], h);
        fcs[c] = std::abs(R::GetArea(xx, fcn[c]));
        fcv[c] = R::GetLineU(fcn[c], fca[c], h) * m.GetVolume(c);
      } else if (fct[c] == Type::regular) {
        fcv[c] = m.GetVolume(c);
      }
    }
  }

  M& m;
  // nodes
  FieldNode<Scal> fnl_; // level-set
  // faces
  FieldFace<Type> fft_; // face type (0: regular, 1: cut, 2: excluded)
  FieldFace<std::vector<Vect>> ffpoly_; // polygon representing f < 0
  FieldFace<Scal> ffs_; // area for which f > 0
  // cells
  FieldCell<Type> fct_; // cell type (0: regular, 1: cut, 2: excluded)
  FieldCell<Vect> fcn_; // unit outer normal (towards excluded domain)
  FieldCell<Scal> fca_; // plane constant
  FieldCell<Scal> fcs_; // area of polygon
  FieldCell<Scal> fcv_; // volume of cut cell
};
