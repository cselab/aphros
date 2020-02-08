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

template <class Scal, size_t N>
static std::array<Scal, N> Mul(
    std::array<Scal, N * N> a, std::array<Scal, N> x) {
  using Int = size_t;
  std::array<Scal, N> r;
  for (Int i = 0; i < N; ++i) {
    r[i] = 0;
    for (Int j = 0; j < N; ++j) {
      r[i] += a[i * N + j] * x[j];
    }
  }
  return r;
}

// Solves linear system a*x=b.
template <class Scal, size_t N>
static std::array<Scal, N> SolveLinear(
    std::array<Scal, N * N> a, std::array<Scal, N> b) {
  using Int = size_t;
  auto aa = [&a](Int i, Int j) -> Scal& { return a[i * N + j]; };
  auto swaprows = [&aa, &b](Int i, Int ip) {
    if (i == ip) {
      return;
    }
    for (Int j = 0; j < N; ++j) {
      std::swap(aa(i, j), aa(ip, j));
    }
    std::swap(b[i], b[ip]);
  };
  auto ipivot = [&aa](const Int j) {
    Int imax = j;
    for (Int i = j + 1; i < N; ++i) {
      if (std::abs(aa(i, j)) > std::abs(aa(imax, j))) {
        imax = i;
      }
    }
    return imax;
  };
  auto addrow = [&aa, &b](Int i, Int ip, Scal ap) {
    for (Int j = 0; j < N; ++j) {
      aa(i, j) += aa(ip, j) * ap;
    }
    b[i] += b[ip] * ap;
  };
  for (Int j = 0; j < N; ++j) {
    const Int ip = ipivot(j);
    swaprows(ip, j);
    for (Int i = j + 1; i < N; ++i) {
      addrow(i, j, -aa(i, j) / aa(j, j));
    }
  }
  std::array<Scal, N> x;
  for (Int i = N; i > 0;) {
    --i;
    Scal t = b[i];
    for (Int j = i + 1; j < N; ++j) {
      t -= aa(i, j) * x[j];
    }
    x[i] = t / aa(i, i);
  }
  return x;
}

// Fits linear function to set of points and values.
template <class Vect, class Scal = typename Vect::value_type>
std::pair<Vect, Scal> FitLinear(
    const std::vector<Vect>& xx, const std::vector<Scal>& uu) {
  assert(xx.size() == uu.size());
  // sum 0.5 * [ (g.dot(x[k]) + u0 - u[k]) ** 2 ] -> min
  using Int = size_t;
  constexpr Int dim = 3;
  constexpr Int N = dim + 1;
  std::array<Scal, N * N> a;
  std::array<Scal, N> b;
  auto aa = [&a](Int i, Int j) -> Scal& { return a[i * N + j]; };
  std::fill(a.begin(), a.end(), 0);
  std::fill(b.begin(), b.end(), 0);
  for (size_t k = 0; k < xx.size(); ++k) {
    for (Int i = 0; i < dim; ++i) {
      for (Int j = 0; j < dim; ++j) {
        aa(i, j) += xx[k][j] * xx[k][i];
      }
      aa(i, dim) += xx[k][i];
      b[i] += uu[k] * xx[k][i];
    }
    for (Int j = 0; j < dim; ++j) {
      aa(dim, j) += xx[k][j];
    }
    aa(dim, dim) += 1;
    b[dim] += uu[k];
  }

  auto x = SolveLinear(a, b);
  return {Vect(x[0], x[1], x[2]), x[3]};
}

template <class Scal>
struct CondEmbed {
  enum class Type { value, gradient };
  Scal u; // value or normal gradient
};

template <class M>
FieldNode<typename M::Scal> InitEmbed(const M& m, const Vars& var, bool verb) {
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
  } else if (name == "list") {
    // TODO revise with bcast
    const std::string fn = var.String["eb_list_path"];
    const size_t edim = var.Int["dim"];
    std::ifstream fin(fn);
    if (verb) {
      std::cout << "Open list of primitives '" << fn << "' for embed"
                << std::endl;
    }
    if (!fin.good()) {
      throw std::runtime_error("Can't open list of primitives");
    }
    auto pp = UPrimList<Scal>::Parse(fin, verb, edim);

    for (auto n : m.AllNodes()) {
      Scal lmax = -std::numeric_limits<Scal>::max(); // maximum level-set
      for (auto& p : pp) {
        lmax = std::max(lmax, p.ls(m.GetNode(n)));
      }
      fnl[n] = lmax;
    }
  } else {
    throw std::runtime_error("Unknown eb_init=" + name);
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
  auto GetFace(IdxCell c, size_t q) const {
    return m.GetFace(c, q);
  }
  auto GetCell(IdxFace f, size_t q) const {
    return m.GetCell(f, q);
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
  // Mid-point interpolation.
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  template <class T>
  FieldEmbed<T> Interpolate(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv) const {
    FieldEmbed<T> feu(m, T(0));
    for (auto f : eb.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const Scal a = 0.5;
      feu[f] = fcu[cm] * (1 - a) + fcu[cp] * a;
    }
    for (auto c : eb.CFaces()) {
      if (bc == 0) {
        feu[c] = bcv;
      } else if (bc == 1) {
        feu[c] = fcu[c] + bcv * GetFaceOffset(c);
      } else {
        throw std::runtime_error(
            "Interpolate: unknown bc=" + std::to_string(bc));
      }
    }
    InterpolateB(fcu, mfc, feu.GetFieldFace(), m);
    return feu;
  }
  // High order gradient-based interpolation with extrapolation
  // boundary conditions.
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  template <class T>
  FieldEmbed<T> Interpolate2(const FieldCell<T>& fcu) const {
    // interpolate from cells to faces with zero-derivative conditions
    auto feu = Interpolate(fcu, MapCondFace(), 1, T(0));

    FieldCell<Vect> fcg(m, Vect(0));

    for (size_t i = 0; i < 20; ++i) {
      // compute gradient from current interpolant
      fcg = Gradient2(feu);

      // goal: compute first-order accurate gradient
      // and using the gradient, do second-order accurate interpolation

      // update feu from fcu and fcg
      for (auto f : eb.Faces()) {
        const IdxCell cm = m.GetCell(f, 0);
        const IdxCell cp = m.GetCell(f, 1);
        if (GetType(cm) == Type::regular && GetType(cp) == Type::regular) {
          feu[f] = (fcu[cm] + fcu[cp]) * 0.5;
        } else {
          const Vect xm = eb.GetCellCenter(cm);
          const Vect xp = eb.GetCellCenter(cp);
          const Vect xf = eb.GetFaceCenter(f);
          feu[f] = (fcu[cm] + fcg[cm].dot(xf - xm) + fcu[cp] +
                    fcg[cp].dot(xf - xp)) *
                   0.5;
        }
      }
      for (auto c : eb.CFaces()) {
        const Vect xc = eb.GetCellCenter(c);
        const Vect xf = eb.GetFaceCenter(c);
        feu[c] = fcu[c] + fcg[c].dot(xf - xc);
      }
    }
    return feu;
  }
  // Upwind interpolation.
  // fcu: field [a]
  // ffv: volume flux
  // mfc: conditions on faces
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  template <class T>
  FieldEmbed<T> InterpolateUpwind(
      const FieldCell<T>& fcu, const FieldEmbed<T>& fev, const MapCondFace& mfc,
      size_t bc, T bcv) const {
    FieldEmbed<T> feu(m, T(0));
    for (auto f : eb.Faces()) {
      const IdxCell c = (fev[f] > 0 ? m.GetCell(f, 0) : m.GetCell(f, 1));
      feu[f] = fcu[c];
    }
    for (auto c : eb.CFaces()) {
      if (fev[c] > 0) {
        feu[c] = fcu[c];
      } else {
        if (bc == 0) {
          feu[c] = bcv;
        } else if (bc == 1) {
          feu[c] = fcu[c] + bcv * GetFaceOffset(c);
        } else {
          throw std::runtime_error(
              "Interpolate: unknown bc=" + std::to_string(bc));
        }
      }
    }
    InterpolateB(fcu, mfc, feu.GetFieldFace(), m);
    return feu;
  }
  // Upwind interpolation.
  // fcu: field [s]
  // fcg: gradient [s]
  // mfc: conditions on faces
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // ffv: volume flux [i]
  // sc: interpolation scheme
  // Returns:
  // field on embedded boundaries [s]
  // TODO: interpolation from upwind cell on embedded boundaries
  //       (see InterpolateUpwind() above)
  FieldEmbed<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const FieldCell<typename M::Vect>& fcg,
      const MapCondFace& mfc, size_t bc, Scal bcv, const FieldFace<Scal>& ffv,
      ConvSc sc) const {
    FieldEmbed<Scal> feu(m, 0.);
    // f = fmm*a[0] + fm*a[1] + fp*a[2]
    const std::array<Scal, 3> a = GetCoeff<Scal>(sc);
    for (auto f : eb.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      if (GetType(cm) == Type::regular && GetType(cp) == Type::regular) {
        if (ffv[f] > 0) {
          feu[f] = 4. * a[0] * fcg[cm].dot(m.GetVectToCell(f, 0)) +
                   a[1] * fcu[cm] + (a[2] + a[0]) * fcu[cp];
        } else if (ffv[f] < 0) {
          feu[f] = 4. * a[0] * fcg[cp].dot(m.GetVectToCell(f, 1)) +
                   a[1] * fcu[cp] + (a[2] + a[0]) * fcu[cm];
        } else {
          feu[f] = (fcu[cm] + fcu[cp]) * 0.5;
        }
      } else {
        const IdxCell c = (ffv[f] > 0 ? cm : cp);
        feu[f] = fcu[c];
      }
    }
    for (auto c : eb.CFaces()) {
      if (bc == 0) {
        feu[c] = bcv;
      } else if (bc == 1) {
        feu[c] = fcu[c] + bcv * GetFaceOffset(c);
      } else {
        throw std::runtime_error(
            "Interpolate: unknown bc=" + std::to_string(bc));
      }
    }
    InterpolateB(fcu, mfc, feu.GetFieldFace(), m);
    return feu;
  }
  Scal ClipGradDenom(Scal dn) const {
    return (dn > 0 ? 1 : -1) *
           std::max(std::abs(dn), m.GetCellSize()[0] * gradlim_);
  }
  // feu: field on embedded boundaries [a]
  // Returns:
  // gradient on cells [a]
  FieldCell<Vect> Gradient(const FieldEmbed<Scal>& feu) const {
    auto& eb = *this;
    FieldCell<Vect> fcg(m, Vect(0));
    for (auto c : m.AllCells()) {
      switch (fct_[c]) {
        case Type::regular: {
          Vect sum(0);
          for (auto q : m.Nci(c)) {
            const IdxFace f = m.GetFace(c, q);
            sum += m.GetOutwardSurface(c, q) * feu[f];
          }
          fcg[c] = sum / m.GetVolume(c);
          break;
        }
        case Type::cut: {
          Vect sum = eb.GetSurface(c) * feu[c];
          for (auto q : eb.Nci(c)) {
            const IdxFace f = m.GetFace(c, q);
            sum += eb.GetOutwardSurface(c, q) * feu[f];
          }
          fcg[c] = sum / eb.GetVolume(c);
          break;
        }
        case Type::excluded:
          break;
      }
    }
    return fcg;
  }
  // feu: field on embedded boundaries [a]
  // Returns:
  // gradient on cells [a]
  FieldCell<Vect> Gradient2(const FieldEmbed<Scal>& feu) const {
    auto& eb = *this;
    FieldCell<Vect> fcg(m, Vect(0));
    for (auto c : m.AllCells()) {
      switch (fct_[c]) {
        case Type::regular: {
          Vect sum(0);
          for (auto q : m.Nci(c)) {
            const IdxFace f = m.GetFace(c, q);
            sum += m.GetOutwardSurface(c, q) * feu[f];
          }
          fcg[c] = sum / m.GetVolume(c);
          break;
        }
        case Type::cut: {
          std::vector<Vect> xx;
          std::vector<Scal> uu;
          xx.push_back(eb.GetFaceCenter(c));
          uu.push_back(feu[c]);
          for (auto q : eb.Nci(c)) {
            const IdxFace f = m.GetFace(c, q);
            xx.push_back(eb.GetFaceCenter(f));
            uu.push_back(feu[f]);
          }
          auto p = FitLinear(xx, uu);
          fcg[c] = p.first;
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
    for (auto c : eb.Cells()) {
      const Scal v0 = m.GetVolume(c);
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
  FieldEmbed<T> Gradient(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv) const {
    FieldEmbed<T> feu(m, T(0)); // FIXME should be nan
    for (auto f : eb.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const Scal dn = ClipGradDenom(
          (eb.GetNormal(f)).dot(eb.GetCellCenter(cp) - eb.GetCellCenter(cm)));
      feu[f] = (fcu[cp] - fcu[cm]) / dn;
    }
    for (auto c : eb.CFaces()) {
      if (bc == 0) {
        const Scal dn = ClipGradDenom(eb.GetFaceOffset(c));
        feu[c] = (bcv - fcu[c]) / dn;
      } else if (bc == 1) {
        feu[c] = bcv;
      } else {
        throw std::runtime_error("Gradient: unknown bc=" + std::to_string(bc));
      }
    }
    GradientB(fcu, mfc, m, feu.GetFieldFace());
    return feu;
  }
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  FieldEmbed<T> Gradient2(const FieldCell<T>& fcu) const {
    auto feu = Interpolate2(fcu);
    auto fcg = Gradient2(feu);
    FieldEmbed<T> feg(m, T(0));
    for (auto f : eb.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      if (GetType(cm) == Type::regular && GetType(cp) == Type::regular) {
        feg[f] = (fcu[cp] - fcu[cm]) / m.GetCellSize()[0];
      } else {
        feg[f] = (fcg[cm] + fcg[cp]).dot(eb.GetNormal(f)) * 0.5;
      }
    }
    for (auto c : eb.CFaces()) {
      feg[c] = fcg[c].dot(eb.GetNormal(c));
    }
    return feg;
  }
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  FieldEmbed<T> Gradient3(const FieldCell<T>& fcu) const {
    FieldEmbed<T> feg(m, T(0));
    for (auto f : eb.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const Scal h = m.GetCellSize()[0];
      if (GetType(cm) == Type::regular && GetType(cp) == Type::regular) {
        feg[f] = (fcu[cp] - fcu[cm]) / h;
      } else {
        const Vect n = eb.GetNormal(cm);
        const size_t qz = m.GetNci(cm, f); // f == GetFace(cm, qz)
        // cp == GetCell(cm, qz)
        //                                             //
        //            ---------------------            //
        //            |         |         |            //
        //            |   cmy   |   cmxy  |            //
        //            |         |         |            //
        //            |\--------|---------|            //
        //            | \       |         |            //
        //            |  \  cm  |   cmx   |            //
        //            |   \     |         |            //
        //            -----\---------------            //
        //                                             //
        const size_t dz = qz / 2; // face direction
        const size_t dx = (dz + 1) % 3; // other directions
        const size_t dy = (dz + 2) % 3; // other directions
        const size_t qx = dx * 2 + (n[dx] < 0 ? 1 : 0);
        const size_t qy = dy * 2 + (n[dy] < 0 ? 1 : 0);
        const IdxCell cmx = m.GetCell(cm, qx);
        const IdxCell cmy = m.GetCell(cm, qy);
        const IdxCell cmxy = m.GetCell(cmx, qy);
        const Scal tx =
            1 - std::abs(eb.GetFaceCenter(f)[dx] - m.GetCenter(f)[dx]) / h;
        const Scal ty =
            1 - std::abs(eb.GetFaceCenter(f)[dy] - m.GetCenter(f)[dy]) / h;
        const Scal g = (fcu[m.GetCell(cm, qz)] - fcu[cm]) / h;
        const Scal gx = (fcu[m.GetCell(cmx, qz)] - fcu[cmx]) / h;
        const Scal gy = (fcu[m.GetCell(cmy, qz)] - fcu[cmy]) / h;
        const Scal gxy = (fcu[m.GetCell(cmxy, qz)] - fcu[cmxy]) / h;
        const Scal gy0 = gx * tx + g * (1 - tx);
        const Scal gy1 = gxy * tx + gy * (1 - tx);
        feg[f] = gy0 * ty + gy1 * (1 - ty);
      }
    }
    for (auto c : eb.CFaces()) {
      feg[c] = 0;
    }
    return feg;
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
