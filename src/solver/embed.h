// Created by Petr Karnakov on 13.05.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <iomanip>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "approx.h"
#include "dump/dumper.h"
#include "dump/vtk.h"
#include "func/primlist.h"
#include "geom/filter.h"
#include "geom/transform.h"
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

  explicit FieldEmbed(const std::pair<GRange<IdxCell>, GRange<IdxFace>>& p)
      : dc_(p.first), df_(p.second) {}
  explicit FieldEmbed(
      const std::pair<GRange<IdxCell>, GRange<IdxFace>>& p, const Value& v)
      : dc_(p.first, v), df_(p.second, v) {}

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

template <class Vect>
const FieldEmbed<typename Vect::value_type> GetComponent(
    const FieldEmbed<Vect>& fev, size_t d) {
  return FieldEmbed<typename Vect::value_type>(
      GetComponent(fev.GetFieldCell(), d), GetComponent(fev.GetFieldFace(), d));
}

template <class T>
class MapEmbed {
 public:
  using Value = T;
  using Container = std::unordered_map<size_t, Value>;

  MapEmbed() = default;
  MapEmbed(const MapEmbed&) = default;
  MapEmbed(MapEmbed&&) = default;
  MapEmbed& operator=(const MapEmbed&) = default;
  MapEmbed& operator=(MapEmbed&&) = default;

  size_t size() const {
    return dc_.size() + df_.size();
  }
  void clear() {
    dc_.clear();
    df_.clear();
  }
  Value& operator[](const IdxCell& c) {
    return dc_[c];
  }
  const Value& at(const IdxCell& c) const {
    return dc_.at(c);
  }
  Value& operator[](const IdxFace& f) {
    return df_[f];
  }
  const Value& at(const IdxFace& f) const {
    return df_.at(f);
  }
  Value* find(IdxCell c) {
    auto it = dc_.find(c);
    return it != dc_.end() ? &it->second : nullptr;
  }
  const Value* find(IdxCell c) const {
    auto it = dc_.find(c);
    return it != dc_.end() ? &it->second : nullptr;
  }
  Value* find(IdxFace f) {
    auto it = df_.find(f);
    return it != df_.end() ? &it->second : nullptr;
  }
  const Value* find(IdxFace f) const {
    auto it = df_.find(f);
    return it != df_.end() ? &it->second : nullptr;
  }
  void erase(const IdxCell& c) {
    dc_.erase(c);
  }
  void erase(const IdxFace& f) {
    df_.erase(f);
  }
  MapCell<Value>& GetMapCell() {
    return dc_;
  }
  const MapCell<Value>& GetMapCell() const {
    return dc_;
  }
  MapFace<Value>& GetMapFace() {
    return df_;
  }
  const MapFace<Value>& GetMapFace() const {
    return df_;
  }

 private:
  MapCell<Value> dc_;
  MapFace<Value> df_;
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
  // Constructor
  // fnl: level-set function on nodes, interface at fnl=0
  Embed(M& m, Scal gradlim) : m(m), eb(*this), gradlim_(gradlim) {}
  Embed(M& m) : Embed(m, 0.5) {}
  template <class Idx>
  operator GRange<Idx>() const {
    return GRange<Idx>(m);
  }
  operator std::pair<GRange<IdxCell>, GRange<IdxFace>>() const {
    return {GRange<IdxCell>(m), GRange<IdxFace>(m)};
  }
  void Init(const FieldNode<Scal>& fnl) {
    auto sem = m.GetSem("init");
    if (sem()) {
      fnl_ = fnl;
      InitFaces(fnl_, fft_, ffpoly_, ffs_, m);
      InitCells(fnl_, ffs_, fct_, fcn_, fca_, fcs_, fcv_, m);
      InitRedistr(fct_, fcv_, fcs_, fc_redistr_, mc_redistr_, m);
      m.Comm(&fc_redistr_);

      // volume of neighbor cells
      fcvst3_.Reinit(m, 0);
      for (auto c : eb.Cells()) {
        for (auto cn : eb.Stencil(c)) {
          fcvst3_[c] += eb.GetVolume(cn);
        }
      }
      m.Comm(&fcvst3_);
    }
  }
  Type GetType(IdxCell c) const {
    return fct_[c];
  }
  Type GetType(IdxFace f) const {
    return fft_[f];
  }
  // Cell indices of cells with embedded boundaries.
  auto CFaces() const {
    return MakeFilterIterator(
        m.Cells(), [this](IdxCell c) { return GetType(c) == Type::cut; });
  }
  auto SuCFaces() const {
    return MakeFilterIterator(
        m.SuCells(), [this](IdxCell c) { return GetType(c) == Type::cut; });
  }
  auto Faces() const {
    return MakeFilterIterator(
        m.Faces(), [this](IdxFace f) { return GetType(f) != Type::excluded; });
  }
  auto SuFaces() const {
    return MakeFilterIterator(m.SuFaces(), [this](IdxFace f) {
      return GetType(f) != Type::excluded;
    });
  }
  IdxCell GetCell(IdxFace f, size_t q) const {
    return m.GetCell(f, q);
  }
  IdxFace GetFace(IdxCell c, size_t q) const {
    return m.GetFace(c, q);
  }
  IdxCell GetCell(IdxCell c, size_t q) const {
    return m.GetCell(c, q);
  }
  auto Nci(IdxCell c) const {
    return MakeFilterIterator(m.Nci(c), [this, c](size_t q) {
      return GetType(m.GetFace(c, q)) != Type::excluded;
    });
  }
  size_t GetNci(IdxCell c, IdxFace f) const {
    return m.GetNci(c, f);
  }
  // Cell indices of non-excluded cells.
  auto Cells() const {
    return MakeFilterIterator(
        m.Cells(), [this](IdxCell c) { return GetType(c) != Type::excluded; });
  }
  // Cell indices of non-excluded cells.
  auto SuCells() const {
    return MakeFilterIterator(m.SuCells(), [this](IdxCell c) {
      return GetType(c) != Type::excluded;
    });
  }
  // Cell indices of non-excluded cells.
  auto AllCells() const {
    return MakeFilterIterator(m.AllCells(), [this](IdxCell c) {
      return GetType(c) != Type::excluded;
    });
  }
  // Cell indices of non-excluded cells from a 3x3x3 stencil.
  auto Stencil(IdxCell c) const {
    return MakeFilterIterator(m.Stencil(c), [this](IdxCell cn) {
      return GetType(cn) != Type::excluded;
    });
  }
  // Cell indices of non-excluded cells from a 3x3x3 stencil.
  auto Stencil5(IdxCell c) const {
    return MakeFilterIterator(m.Stencil5(c), [this](IdxCell cn) {
      return GetType(cn) != Type::excluded;
    });
  }
  M& GetMesh() {
    return m;
  }
  const M& GetMesh() const {
    return m;
  }
  Vect GetCellSize() const {
    return m.GetCellSize();
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
  Scal GetAreaFraction(IdxFace f) const {
    return ffs_[f] / m.GetArea(f);
  }
  Scal GetArea(IdxCell c) const {
    return fcs_[c];
  }
  Vect GetSurface(IdxFace f) const {
    return GetNormal(f) * GetArea(f);
  }
  Vect GetSurface(IdxCell c) const {
    return GetNormal(c) * GetArea(c);
  }
  Vect GetOutwardSurface(IdxCell c, size_t q) const {
    return GetSurface(m.GetFace(c, q)) * m.GetOutwardFactor(c, q);
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
  Scal GetVolumeFraction(IdxCell c) const {
    return GetVolume(c) / m.GetVolume(c);
  }
  Scal GetVolumeStencilSum(IdxCell c) const {
    return fcvst3_[c];
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
  // TODO: cache
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
  std::vector<Vect> GetCutPoly(IdxCell c) const {
    return R::GetCutPoly(m.GetCenter(c), fcn_[c], fca_[c], m.GetCellSize());
  }
  Vect GetNode(IdxNode n) const {
    return m.GetNode(n);
  }
  // TODO: cache
  Vect GetCellCenter(IdxCell c) const {
    switch (fct_[c]) {
      case Type::regular:
        return m.GetCenter(c);
      case Type::cut: {
        const Scal w = GetArea(c);
        Vect sum = GetFaceCenter(c) * w;
        Scal sumw = w;
        for (auto q : this->Nci(c)) {
          auto f = m.GetFace(c, q);
          const Scal w = GetArea(f);
          sum += GetFaceCenter(f) * w;
          sumw += w;
        }
        return sum / sumw;
      }
      default:
        return GetNan<Vect>();
    }
  }
  Scal GetRedistr(IdxCell c) const {
    return fc_redistr_[c];
  }
  const FieldCell<Scal>& GetRedistr() const {
    return fc_redistr_;
  }
  const std::vector<std::pair<IdxCell, Scal>>& GetRedistrList(IdxCell c) const {
    return *mc_redistr_.find(c);
  }
  Scal ClipGradDenom(Scal dn) const {
    return (dn > 0 ? 1 : -1) *
           std::max(std::abs(dn), m.GetCellSize()[0] * gradlim_);
  }

  void DumpPoly(std::string filename, bool vtkbin, bool vtkmerge) const {
    DumpPoly(
        filename, ffs_, fft_, fcs_, fct_, fcn_, fca_, ffpoly_, vtkbin, vtkmerge,
        m);
  }
  void DumpPoly() const {
    DumpPoly("eb.vtk", true, true);
  }
  void DumpPoly(bool vtkbin, bool vtkmerge) const {
    DumpPoly("eb.vtk", vtkbin, vtkmerge);
  }
  // xc: plane origin
  // n: plane normal
  void DumpPlaneSection(std::string filename, Vect xc, Vect n) const {
    DumpPlaneSection(filename, fct_, fcn_, fca_, xc, n, m);
  }
  // xc: plane origin
  // n: plane normal
  void DumpPlaneSection(Vect xc, Vect n) const {
    DumpPlaneSection("eb.dat", xc, n);
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
      bool vtkbin, bool vtkmerge, M& m) {
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
      m.Reduce(std::make_shared<TS>(&dls));
      m.Reduce(std::make_shared<TS>(&dlf));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::cout << std::fixed << std::setprecision(8) << "dump"
                  << " to " << fn << std::endl;
        WriteVtkPoly<Vect>(
            fn, dl, nullptr, {&dld, &dls, &dlf}, {"dir", "area", "face"},
            "Embedded boundary", true, vtkbin, vtkmerge);
      }
    }
  }
  // Writes line segments of cross-section of cut cells
  // fct: cell types
  // fcn: normals
  // fca: plane constants
  // plane_x: plane origin
  // plane_n: plane normal
  static void DumpPlaneSection(
      const std::string fn, const FieldCell<Type>& fct,
      const FieldCell<Vect>& fcn, const FieldCell<Scal>& fca, Vect plane_x,
      Vect plane_n, M& m) {
    auto sem = m.GetSem("dumpplanesection");
    struct {
      std::vector<std::vector<Vect>> dl;
    } * ctx(sem);
    auto& dl = ctx->dl;
    if (sem("local")) {
      for (auto c : m.Cells()) {
        if (fct[c] == Type::cut) {
          auto xx =
              R::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], m.GetCellSize());
          std::array<Vect, 2> e; // ends of intersection
          if (R::GetInterPoly(xx, plane_x, plane_n, e)) {
            dl.push_back({e[0], e[1]});
          }
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::cout << std::fixed << std::setprecision(8) << "dump"
                  << " to " << fn << std::endl;
        std::ofstream f(fn);
        for (auto e : dl) {
          f << e[0][0] << ' ' << e[0][1] << ' ' << e[0][2] << ' ';
          f << e[1][0] << ' ' << e[1][1] << ' ' << e[1][2] << '\n';
        }
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
  // Area of plane convex polygon.
  // xx: points of polygon
  // n: normal to plane
  // Returns:
  // a: area, positive if <xx[0]-xc,xx[1]-xc,n> is positively oriented
  static Scal GetArea0(const std::vector<Vect>& xx, Vect n) {
    size_t sx = xx.size();

    if (!sx) {
      return 0.;
    }
    // unit normal
    n /= n.norm();

    Scal a = 0.;
    for (size_t i = 1; i < sx; ++i) {
      const size_t ip = (i + 1 == sx ? 0 : i + 1);
      a += (xx[i] - xx[0]).cross(xx[ip] - xx[0]).dot(n) * 0.5;
    }
    return a;
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
        const size_t ep = (e + 1) % em;
        const IdxNode n = m.GetNode(f, e);
        const IdxNode np = m.GetNode(f, ep);
        const Scal l = fnl[n];
        const Scal lp = fnl[np];
        const Vect x = m.GetNode(n);
        const Vect xp = m.GetNode(np);
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
        Vect sum(0);
        for (auto q : m.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          sum += m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q);
        }
        fcs[c] = -sum.dot(fcn[c]);

        fcv[c] = R::GetLineU(fcn[c], fca[c], h) * m.GetVolume(c);
      } else if (fct[c] == Type::regular) {
        fcv[c] = m.GetVolume(c);
      }
    }
  }
  static void InitRedistr(
      const FieldCell<Type>& fct, const FieldCell<Scal>& fcv,
      const FieldCell<Scal>& fcs, FieldCell<Scal>& fc_redistr,
      MapCell<std::vector<std::pair<IdxCell, Scal>>>& mc_redistr, const M& m) {
    fc_redistr.Reinit(m, 1);
    const Scal hreg = m.GetCellSize()[0]; // XXX adhoc cubic cells
    for (auto c : m.Cells()) {
      if (fct[c] == Type::cut) {
        const Scal a = std::min(1., fcv[c] / (hreg * fcs[c]));
        std::vector<std::pair<IdxCell, Scal>> vp;
        for (auto cn : m.Stencil(c)) {
          if (fct[cn] != Type::excluded && cn != c) {
            vp.push_back({cn, 1});
          }
        }
        Scal sum = 0;
        for (auto& p : vp) {
          sum += p.second * fcv[p.first];
        }
        for (auto& p : vp) {
          p.second *= fcv[c] * (1 - a) / sum;
        }
        fc_redistr[c] = a;
        mc_redistr[c] = vp;
      }
    }
  }

  M& m;
  const Embed& eb;
  const Scal gradlim_;
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
  // redistribution coefficients
  MapCell<std::vector<std::pair<IdxCell, Scal>>> mc_redistr_;
  FieldCell<Scal> fc_redistr_;
  // volume of neighbor cells
  FieldCell<Scal> fcvst3_; // volume of neighbors in stencil 3x3x3
};
