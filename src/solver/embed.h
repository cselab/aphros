// Created by Petr Karnakov on 13.05.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <iomanip>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "geom/filter.h"
#include "geom/mesh.h"
#include "reconst.h"

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
  Value& at(const IdxCell& c) {
    return dc_.at(c);
  }
  const Value& at(const IdxCell& c) const {
    return dc_.at(c);
  }
  Value& operator[](const IdxFace& f) {
    return df_[f];
  }
  Value& at(const IdxFace& f) {
    return df_.at(f);
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
  template <class F>
  void LoopPairs(F lambda) const {
    for (auto& p : df_) {
      lambda(p);
    }
    for (auto& p : dc_) {
      lambda(p);
    }
  }
  // Iterates over faces and passes the index of neighbor cell.
  // Value should have member variable `nci`.
  // F: function void(IdxFace or IdxCell, IdxCell, BC)
  template <class EB, class F>
  void LoopBCond(const EB& eb, F lambda) {
    for (auto& p : df_) {
      const IdxFace f = p.first;
      const auto& bc = p.second;
      lambda(f, eb.GetCell(f, bc.nci), bc);
    }
    for (auto& p : dc_) {
      const IdxCell c = p.first;
      const auto& bc = p.second;
      lambda(c, c, bc);
    }
  }
  template <class EB, class F>
  void LoopBCond(const EB& eb, F lambda) const {
    for (auto& p : df_) {
      const IdxFace f = p.first;
      const auto& bc = p.second;
      lambda(f, eb.GetCell(f, bc.nci), bc);
    }
    for (auto& p : dc_) {
      const IdxCell c = p.first;
      const auto& bc = p.second;
      lambda(c, c, bc);
    }
  }

 private:
  MapCell<Value> dc_;
  MapFace<Value> df_;
};

// Embedded boundaries.
template <class M_>
class Embed {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

 private:
  using R = Reconst<Scal>;

 public:

  enum class Type { regular, cut, excluded };
  class NciEmbed {};
  // Constructor
  // fnl: level-set function on nodes, interface at fnl=0
  Embed(M& m, Scal gradlim) : m(m), eb(*this), gradlim_(gradlim) {}
  explicit Embed(M& m) : Embed(m, 0.5) {}
  // Initializes embedded boundaries with level-set function fnl.
  // Suspendable, requires communication.
  void Init(const FieldNode<Scal>& fnl);
  template <class Idx>
  operator GRange<Idx>() const {
    return GRange<Idx>(m);
  }
  operator std::pair<GRange<IdxCell>, GRange<IdxFace>>() const {
    return {GRange<IdxCell>(m), GRange<IdxFace>(m)};
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
  template <class F>
  void LoopFaces(F lambda) const {
    for (auto c : CFaces()) {
      lambda(c);
    }
    for (auto f : Faces()) {
      lambda(f);
    }
  }
  template <class F>
  void LoopSuFaces(F lambda) const {
    for (auto c : SuCFaces()) {
      lambda(c);
    }
    for (auto f : SuFaces()) {
      lambda(f);
    }
  }
  template <class F>
  void LoopNci(IdxCell c, F lambda) const {
    for (auto q : Nci(c)) {
      lambda(q);
    }
    lambda(NciEmbed());
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
  IdxCell GetFace(IdxCell c, NciEmbed) const {
    return c;
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
  Scal GetOutwardFactor(IdxCell c, size_t q) const {
    return m.GetOutwardFactor(c, q);
  }
  Scal GetOutwardFactor(IdxCell, NciEmbed) const {
    return 1;
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
  // TODO: cache
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
      FieldFace<std::vector<Vect>>& ffpoly, FieldFace<Scal>& ffs, const M& m);
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
      FieldCell<Scal>& fcs, FieldCell<Scal>& fcv, const M& m);
  static void InitRedistr(
      const FieldCell<Type>& fct, const FieldCell<Scal>& fcv,
      const FieldCell<Scal>& fcs, FieldCell<Scal>& fc_redistr,
      MapCell<std::vector<std::pair<IdxCell, Scal>>>& mc_redistr, const M& m);
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
      bool vtkbin, bool vtkmerge, M& m);
  // Writes line segments of cross-section of cut cells
  // fct: cell types
  // fcn: normals
  // fca: plane constants
  // plane_x: plane origin
  // plane_n: plane normal
  static void DumpPlaneSection(
      const std::string fn, const FieldCell<Type>& fct,
      const FieldCell<Vect>& fcn, const FieldCell<Scal>& fca, Vect plane_x,
      Vect plane_n, M& m);

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

template <class M>
struct EmbedTraits {
  template <class T>
  using FieldFaceb = FieldFace<T>;
};

template <class M>
struct EmbedTraits<Embed<M>> {
  template <class T>
  using FieldFaceb = FieldEmbed<T>;
};
