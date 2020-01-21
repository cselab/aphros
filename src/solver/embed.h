// Created by Petr Karnakov on 13.05.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <iomanip>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "dump/dumper.h"
#include "dump/vtk.h"
#include "geom/filter.h"
#include "geom/transform.h"
#include "reconst.h"
#include "solver.h"

template <class M>
FieldNode<typename M::Scal> InitEmbed(const M& m, const Vars& var) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  FieldNode<Scal> fnl(m); // level-set
  const auto name = var.String["eb_init"];
  if (name == "none") {
    fnl.Reinit(m, 1);
  } else if (name == "box") {
    const Vect xc(var.Vect["eb_box_c"]);
    const Vect r(var.Vect["eb_box_r"]);
    const Scal angle = M_PI * var.Double["eb_box_angle"];
    for (auto n : m.AllNodes()) {
      const Vect x = m.GetNode(n);
      auto rot = [angle](Vect xx) {
        const Scal sin = std::sin(angle);
        const Scal cos = std::cos(angle);
        const Scal x = xx[0];
        const Scal y = xx[1];
        const Scal z = xx[2];
        return Vect(x * cos - y * sin, x * sin + y * cos, z);
      };
      fnl[n] = (1 - (rot(x - xc) / r).norminf()) * (r / m.GetCellSize()).min();
    }
  } else if (name == "sphere") {
    const Vect xc(var.Vect["eb_sphere_c"]);
    const Vect r(var.Vect["eb_sphere_r"]);
    const Scal angle = M_PI * var.Double["eb_sphere_angle"];
    for (auto n : m.AllNodes()) {
      const Vect x = m.GetNode(n);
      auto rot = [angle](Vect xx) {
        const Scal sin = std::sin(angle);
        const Scal cos = std::cos(angle);
        const Scal x = xx[0];
        const Scal y = xx[1];
        const Scal z = xx[2];
        return Vect(x * cos - y * sin, x * sin + y * cos, z);
      };
      fnl[n] = (rot(x - xc) / r).norm() - 1;
    }
  }
  if (var.Int["eb_init_inverse"]) {
    for (auto n : m.AllNodes()) {
      fnl[n] = -fnl[n];
    }
  }
  return fnl;
}

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

/*
// Idx_: instance of GIdx
template <class T>
class MapEmbed {
 public:
  using Idx = Idx_;
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
    return dc_[c.GetRaw()];
  }
  const Value& operator[](const IdxCell& c) const {
    return dc_[c.GetRaw()];
  }
  Value& operator[](const IdxFace& f) {
    return df_[f.GetRaw()];
  }
  const Value& operator[](const IdxFace& f) const {
    return df_[f.GetRaw()];
  }
  Value* find(IdxCell c) {
    auto it = dc_.find(c.GetRaw());
    return it != dc_.end() ? &it->second : nullptr;
  }
  const Value* find(IdxCell c) const {
    auto it = dc_.find(c.GetRaw());
    return it != dc_.end() ? &it->second : nullptr;
  }
  Value* find(IdxFace f) {
    auto it = df_.find(f.GetRaw());
    return it != df_.end() ? &it->second : nullptr;
  }
  const Value* find(IdxFace f) const {
    auto it = df_.find(f.GetRaw());
    return it != df_.end() ? &it->second : nullptr;
  }
  void erase(const IdxCell& c) {
    dc_.erase(c.GetRaw());
  }
  void erase(const IdxFace& f) {
    df_.erase(f.GetRaw());
  }

 private:
  Cont d_;
};
*/

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
  Embed(M& m, const FieldNode<Scal>& fnl) : m(m), eb(*this), fnl_(fnl) {
    InitFaces(fnl_, fft_, ffpoly_, ffs_, m);
    InitCells(fnl_, ffs_, fct_, fcn_, fca_, fcs_, fcv_, m);
    InitRedistr(fct_, fcv_, fcs_, fc_redistr_, mc_redistr_, m);

    fcvst3_.Reinit(m, 0);
    for (auto c : eb.Cells()) {
      for (auto cn : eb.Stencil(c)) {
        fcvst3_[c] += eb.GetVolume(cn);
      }
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
  auto Faces() const {
    return MakeFilterIterator(
        m.Faces(), [this](IdxFace f) { return GetType(f) != Type::excluded; });
  }
  auto Nci(IdxCell c) const {
    return MakeFilterIterator(m.Nci(c), [this, c](size_t q) {
      return GetType(m.GetFace(c, q)) != Type::excluded;
    });
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
  Vect GetSurface(IdxFace f) const {
    return GetNormal(f) * GetArea(f);
  }
  Vect GetSurface(IdxCell c) const {
    return GetNormal(c) * GetArea(c);
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

  // feu: field on embedded boundaries [a]
  // Returns:
  // field on cells [a]
  template <class T>
  FieldCell<T> Interpolate(const FieldEmbed<T>& feu) const {
    FieldCell<T> fcu(m, T(0)); // FIXME should be nan
    for (auto c : m.AllCells()) {
      switch (fct_[c]) {
        case Type::regular: {
          T sum(0);
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
          T sum = feu[c] * w;
          Scal sumw = w;
          for (auto q : this->Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            const Scal w = 1 / std::abs(GetFaceOffset(c, q));
            sum += feu[f] * w;
            sumw += w;
          }
          fcu[c] = sum / sumw;
          break;
        }
        case Type::excluded:
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
  template <class T>
  FieldEmbed<T> Interpolate(const FieldCell<T>& fcu, size_t bc, T bcv) const {
    FieldEmbed<T> feu(m, T(0)); // FIXME should be nan
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
            feu[c] = fcu[c] + bcv * GetFaceOffset(c);
          } else {
            throw std::runtime_error(
                "Interpolate: unknown bc=" + std::to_string(bc));
          }
          break;
        }
        case Type::regular:
        case Type::excluded: {
          break;
        }
      }
    }
    return feu;
  }
  // feu: field on embedded boundaries [a]
  // Returns:
  // gradient on cells [a]
  FieldCell<Vect> Gradient(const FieldEmbed<Scal>& feu) const {
    FieldCell<Vect> fcg(m, Vect(0)); // FIXME should be nan
    for (auto c : m.AllCells()) {
      switch (fct_[c]) {
        case Type::regular: {
          Vect sum(0);
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            sum += m.GetOutwardSurface(c, q) * feu[f];
          }
          fcg[c] = sum / m.GetVolume(c);
          break;
        }
        case Type::cut: {
          Vect sum = GetNormal(c) * (GetArea(c) * feu[c]);
          for (auto q : this->Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            sum +=
                GetNormal(f) * (m.GetOutwardFactor(c, q) * GetArea(f) * feu[f]);
          }
          fcg[c] = sum / GetVolume(c);
          break;
        }
        case Type::excluded:
          break;
      }
    }
    return fcg;
  }
  template <class T>
  FieldCell<T> AverageCutCells(const FieldCell<T>& fcu) const {
    FieldCell<T> fcr = fcu;
    for (auto c : eb.Cells()) {
      if (eb.GetType(c) == Type::cut) {
        const Scal v = eb.GetVolume(c);
        T sum = fcu[c] * v;
        Scal sumv = v;
        for (IdxCell cn : eb.Stencil(c)) {
          const Scal vn = eb.GetVolume(cn);
          sum += fcu[cn] * vn;
          sumv += vn;
        }
        fcr[c] = sum / sumv;
      }
    }
    return fcr;
  }
  template <class T>
  FieldCell<T> RedistributeCutCells(const FieldCell<T>& fcu) const {
    FieldCell<T> fcr = fcu;
    const Scal v0 = m.GetCellSize().prod();
    for (auto c : eb.Cells()) {
      const Scal v = eb.GetVolume(c);
      // excess quantity
      const T du = fcu[c] * (1 - v / v0);
      // subtract from current cell
      fcr[c] -= du;
      // add from neighbor cells proportional to their volume
      for (auto cn : eb.Stencil(c)) {
        if (c != cn) {
          const Scal vn = eb.GetVolume(cn);
          // excess quantity in cell cn
          const T dun = fcu[cn] * (1 - vn / v0);
          fcr[c] += dun * (v / (fcvst3_[cn] - vn));
        }
      }
    }
    return fcr;
  }
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  FieldEmbed<T> Gradient(const FieldCell<T>& fcu, size_t bc, T bcv) const {
    FieldEmbed<T> feu(m, T(0)); // FIXME should be nan
    for (auto f : eb.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      Scal dn =
          (eb.GetCellCenter(cp) - eb.GetCellCenter(cm)).dot(eb.GetNormal(f));
      dn = (dn > 0 ? 1 : -1) * std::max(std::abs(dn), m.GetCellSize()[0] * 0.5);
      feu[f] = (fcu[cp] - fcu[cm]) / dn;
    }
    for (auto c : eb.CFaces()) {
      if (bc == 0) {
        const Scal dn = std::max(eb.GetFaceOffset(c), m.GetCellSize()[0] * 0.5);
        feu[c] = (bcv - fcu[c]) / dn;
      } else if (bc == 1) {
        feu[c] = bcv;
      } else {
        throw std::runtime_error("Gradient: unknown bc=" + std::to_string(bc));
      }
    }
    return feu;
  }
  void DumpPoly(std::string filename, bool vtkbin, bool vtkmerge) const {
    DumpPoly(
        filename, ffs_, fft_, fcs_, fct_, fcn_, fca_, ffpoly_, vtkbin, vtkmerge,
        m);
  }
  void DumpPoly() const {
    DumpPoly(
        "eb.vtk", ffs_, fft_, fcs_, fct_, fcn_, fca_, ffpoly_, true, true, m);
  }
  void DumpPoly(bool vtkbin, bool vtkmerge) const {
    DumpPoly(
        "eb.vtk", ffs_, fft_, fcs_, fct_, fcn_, fca_, ffpoly_, vtkbin, vtkmerge,
        m);
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
  static void InitRedistr(
      const FieldCell<Type>& fct, const FieldCell<Scal>& fcv,
      const FieldCell<Scal>& fcs, FieldCell<Scal>& fc_redistr,
      MapCell<std::vector<std::pair<IdxCell, Scal>>>& mc_redistr, const M& m) {
    fc_redistr.Reinit(m, 1);
    const int sw = 1; // stencil halfwidth
    using MIdx = typename M::MIdx;
    auto& bc = m.GetIndexCells();
    GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sw * 2 + 1)); // offsets
    const Scal hreg = m.GetCellSize()[0]; // XXX adhoc cubic cells
    for (auto c : m.Cells()) {
      if (fct[c] == Type::cut) {
        const Scal a = std::min(1., fcv[c] / (hreg * fcs[c]));
        std::vector<std::pair<IdxCell, Scal>> vp;
        const MIdx w = bc.GetMIdx(c);
        for (auto wo : bo) {
          const MIdx wn = w + wo;
          const IdxCell cn = bc.GetIdx(wn);
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
