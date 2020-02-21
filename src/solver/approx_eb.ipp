// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "approx_eb.h"

#include "approx.h"
#include "func/primlist.h"

template <class Scal>
auto ULinear<Scal>::FitLinear(
    const std::vector<Vect>& xx, const std::vector<Scal>& uu)
    -> std::pair<Vect, Scal> {
  assert(xx.size() == uu.size());
  // sum 0.5 * [ (g.dot(x[k]) + u0 - u[k]) ** 2 ] -> min
  using Int = size_t;
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

  auto v = SolveLinear(a, b);
  return {Vect(v[0], v[1], v[2]), v[3]};
}

template <class Scal>
auto ULinear<Scal>::FitLinear(
    const std::vector<Vect>& xx, const std::vector<Vect>& uu)
    -> std::pair<generic::Vect<Vect, dim>, Vect> {
  std::pair<generic::Vect<Vect, dim>, Vect> p;
  for (size_t d = 0; d < dim; ++d) {
    std::vector<Scal> uud;
    for (auto u : uu) {
      uud.push_back(u[d]);
    }
    auto pd = FitLinear(xx, uud);
    p.second[d] = pd.second;
    p.first[0][d] = pd.first[0];
    p.first[1][d] = pd.first[1];
    p.first[2][d] = pd.first[2];
  }
  return p;
}

template <class M>
auto UEmbed<M>::InitEmbed(const M& m, const Vars& var, bool verb)
    -> FieldNode<Scal> {
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

    auto lsmax = [&pp](Vect x) {
      Scal lmax = -std::numeric_limits<Scal>::max(); // maximum level-set
      for (size_t i = 0; i < pp.size(); ++i) {
        const auto& p = pp[i];
        Scal li = p.ls(x);
        if (p.mod_minus) {
          li = -li;
        }
        if (p.mod_and) {
          lmax = std::min(lmax, li);
        } else {
          if (li > lmax) {
            lmax = li;
          }
        }
      }
      return lmax;
    };
    for (auto n : m.AllNodes()) {
      fnl[n] = lsmax(m.GetNode(n));
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

template <class M, class T>
auto FitLinear(IdxCell c, const FieldCell<T>& fcu, const Embed<M>& eb) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto& m = eb.GetMesh();
  std::vector<Vect> xx;
  std::vector<T> uu;
  for (auto cn : eb.Stencil(c)) {
    xx.push_back(m.GetCenter(cn));
    uu.push_back(fcu[cn]);
  }
  return ULinear<Scal>::FitLinear(xx, uu);
}

template <class M, class T>
auto EvalLinearFit(
    typename M::Vect x, IdxCell c, const FieldCell<T>& fcu,
    const Embed<M>& eb) {
  auto p = FitLinear(c, fcu, eb);
  return ULinear<typename M::Scal>::EvalLinear(p, x);
}

template <class Scal>
struct CombineMixed {
  using Vect = generic::Vect<Scal, 3>;
  Scal operator()(Scal un, Scal, Vect) {
    return un;
  }
  // un: normal component
  // ut: tangential component
  // n: normal
  Vect operator()(Vect un, Vect ut, Vect n) {
    return un.proj(n) + ut.orth(n);
  }
};

template <class M>
template <class T>
auto UEmbed<M>::Interpolate(
    const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const EB& eb)
    -> FieldEmbed<T> {
  FieldEmbed<T> feu(eb, T(0));

  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    feu[f] = (fcu[cp] + fcu[cm]) * 0.5;
  }
  feu.GetFieldFace() = InterpolateBilinearFaces(feu.GetFieldFace(), eb);

  auto calc = [&](auto i, IdxCell c, auto& bc) {
    const Scal h = eb.GetCellSize()[0];
    const T& val = bc.val;
    switch (bc.type) {
      case BCondType::dirichlet:
        return val;
      case BCondType::neumann: {
        const Vect x = eb.GetFaceCenter(i) - eb.GetNormal(i) * h;
        const T u = EvalLinearFit(x, c, fcu, eb);
        return u + val * h;
      }
      case BCondType::mixed: {
        const Vect x = eb.GetFaceCenter(i) - eb.GetNormal(i) * h;
        const T u = EvalLinearFit(x, c, fcu, eb);
        return CombineMixed<Scal>()(val, u + val * h, eb.GetNormal(c));
      }
      case BCondType::extrap:
        return EvalLinearFit(eb.GetFaceCenter(i), c, fcu, eb);
    }
    return T(0);
  };
  // embedded faces
  for (auto& p : mebc.GetMapCell()) {
    const IdxCell c = p.first;
    const auto& bc = p.second;
    feu[c] = calc(c, c, bc);
  }
  // faces with boundary conditions
  for (auto& p : mebc.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    feu[f] = calc(f, eb.GetCell(f, bc.nci), bc);
  }
  return feu;
}

template <class M>
template <class T>
auto UEmbed<M>::Interpolate(const FieldEmbed<T>& feu, const EB& eb)
    -> FieldCell<T> {
  auto& m = eb.GetMesh();
  FieldCell<T> fcu(eb, T(0));
  for (auto c : eb.AllCells()) {
    if (eb.GetType(c) == Type::regular) {
      T sum(0);
      Scal sumw = 0;
      for (auto q : m.Nci(c)) {
        IdxFace f = eb.GetFace(c, q);
        sum += feu[f];
        sumw += 1.;
      }
      fcu[c] = sum / sumw;
    } else { // Type::cut
      const Scal w = 1 / std::abs(eb.GetFaceOffset(c));
      T sum = feu[c] * w;
      Scal sumw = w;
      for (auto q : eb.Nci(c)) {
        IdxFace f = eb.GetFace(c, q);
        const Scal w = 1 / std::abs(eb.GetFaceOffset(c, q));
        sum += feu[f] * w;
        sumw += w;
      }
      fcu[c] = sum / sumw;
    }
  }
  return fcu;
}

template <class M>
template <class T>
auto UEmbed<M>::Interpolate(
    const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
    const EB& eb) -> FieldEmbed<T> {
  FieldEmbed<T> feu(eb, T(0));
  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    const Scal a = 0.5;
    feu[f] = fcu[cm] * (1 - a) + fcu[cp] * a;
  }
  for (auto c : eb.SuCFaces()) {
    if (bc == 0) {
      feu[c] = bcv;
    } else if (bc == 1) {
      feu[c] = fcu[c] + bcv * eb.GetFaceOffset(c);
    } else {
      throw std::runtime_error("Interpolate: unknown bc=" + std::to_string(bc));
    }
  }
  InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
  return feu;
}

template <class M>
template <class T>
auto UEmbed<M>::InterpolateUpwind(
    const FieldCell<T>& fcu, const FieldEmbed<Scal>& fev,
    const MapCondFace& mfc, size_t bc, T bcv, const EB& eb) -> FieldEmbed<T> {
  FieldEmbed<T> feu(eb, T(0));
  for (auto f : eb.SuFaces()) {
    const IdxCell c = (fev[f] > 0 ? eb.GetCell(f, 0) : eb.GetCell(f, 1));
    feu[f] = fcu[c];
  }
  for (auto c : eb.SuCFaces()) {
    if (fev[c] > 0) {
      feu[c] = fcu[c];
    } else {
      if (bc == 0) {
        feu[c] = bcv;
      } else if (bc == 1) {
        feu[c] = fcu[c] + bcv * eb.GetFaceOffset(c);
      } else {
        throw std::runtime_error(
            "Interpolate: unknown bc=" + std::to_string(bc));
      }
    }
  }
  InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
  return feu;
}

template <class M>
auto UEmbed<M>::InterpolateUpwindBilinear(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
    const FieldFace<Scal>& ffv, ConvSc sc, const EB& eb) -> FieldFace<Scal> {
  auto& m = eb.GetMesh();
  FieldFace<Scal> ffu(eb, 0.);
  // f = fmm*a[0] + fm*a[1] + fp*a[2]
  const std::array<Scal, 3> a = GetCoeff<Scal>(sc);
  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    if (ffv[f] > 0) {
      ffu[f] = 4. * a[0] * fcg[cm].dot(m.GetVectToCell(f, 0)) + a[1] * fcu[cm] +
               (a[2] + a[0]) * fcu[cp];
    } else if (ffv[f] < 0) {
      ffu[f] = 4. * a[0] * fcg[cp].dot(m.GetVectToCell(f, 1)) + a[1] * fcu[cp] +
               (a[2] + a[0]) * fcu[cm];
    } else {
      ffu[f] = (fcu[cm] + fcu[cp]) * 0.5;
    }
  }
  return InterpolateBilinearFaces(ffu, eb);
}

template <class M>
auto UEmbed<M>::InterpolateUpwindBilinear(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
    const MapCondFace& mfc, size_t bc, const MapCell<Scal>& mcu,
    const FieldEmbed<Scal>& fev, ConvSc sc, const EB& eb) -> FieldEmbed<Scal> {
  FieldEmbed<Scal> feu(eb, 0);
  feu.GetFieldFace() =
      InterpolateUpwindBilinear(fcu, fcg, fev.GetFieldFace(), sc, eb);
  InterpolateUpwindEmbedFaces(fcu, bc, mcu, fev, feu, eb);
  InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
  return feu;
}

template <class M>
auto UEmbed<M>::InterpolateUpwindBilinear(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
    const MapCondFace& mfc, size_t bc, Scal bcv, const FieldEmbed<Scal>& fev,
    ConvSc sc, const EB& eb) -> FieldEmbed<Scal> {
  MapCell<Scal> mcu;
  for (auto c : eb.SuCFaces()) {
    mcu[c] = bcv;
  }
  return InterpolateUpwindBilinear(fcu, fcg, mfc, bc, mcu, fev, sc, eb);
}

template <class M>
auto UEmbed<M>::InterpolateUpwind(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
    const MapCondFace& mfc, size_t bc, const MapCell<Scal>& mcu,
    const FieldEmbed<Scal>& fev, ConvSc sc, const EB& eb) -> FieldEmbed<Scal> {
  auto& m = eb.GetMesh();
  FieldEmbed<Scal> feu(eb, 0.);
  // f = fmm*a[0] + fm*a[1] + fp*a[2]
  const std::array<Scal, 3> a = GetCoeff<Scal>(sc);
  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    if (eb.GetType(cm) == Type::regular && eb.GetType(cp) == Type::regular) {
      if (fev[f] > 0) {
        feu[f] = 4. * a[0] * fcg[cm].dot(m.GetVectToCell(f, 0)) +
                 a[1] * fcu[cm] + (a[2] + a[0]) * fcu[cp];
      } else if (fev[f] < 0) {
        feu[f] = 4. * a[0] * fcg[cp].dot(m.GetVectToCell(f, 1)) +
                 a[1] * fcu[cp] + (a[2] + a[0]) * fcu[cm];
      } else {
        feu[f] = (fcu[cm] + fcu[cp]) * 0.5;
      }
    } else {
      const IdxCell c = (fev[f] > 0 ? cm : cp);
      feu[f] = fcu[c];
    }
  }
  InterpolateUpwindEmbedFaces(fcu, bc, mcu, fev, feu, eb);
  InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
  return feu;
}

template <class M>
auto UEmbed<M>::InterpolateUpwind(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
    const MapCondFace& mfc, size_t bc, Scal bcv, const FieldEmbed<Scal>& fev,
    ConvSc sc, const EB& eb) -> FieldEmbed<Scal> {
  MapCell<Scal> mcu;
  for (auto c : eb.SuCFaces()) {
    mcu[c] = bcv;
  }
  return InterpolateUpwind(fcu, fcg, mfc, bc, mcu, fev, sc, eb);
}

template <class M>
auto UEmbed<M>::Gradient(const FieldEmbed<Scal>& feu, const EB& eb)
    -> FieldCell<Vect> {
  auto& m = eb.GetMesh();
  FieldCell<Vect> fcg(eb, Vect(0));
  for (auto c : eb.AllCells()) {
    if (eb.GetType(c) == Type::regular) {
      Vect sum(0);
      for (auto q : m.Nci(c)) {
        const IdxFace f = eb.GetFace(c, q);
        sum += m.GetOutwardSurface(c, q) * feu[f];
      }
      fcg[c] = sum / m.GetVolume(c);
    } else { // Type::cut
      Vect sum = eb.GetSurface(c) * feu[c];
      for (auto q : eb.Nci(c)) {
        const IdxFace f = eb.GetFace(c, q);
        sum += eb.GetOutwardSurface(c, q) * feu[f];
      }
      fcg[c] = sum / eb.GetVolume(c);
    }
  }
  return fcg;
}

template <class M>
auto UEmbed<M>::GradientLinearFit(const FieldEmbed<Scal>& feu, const EB& eb)
    -> FieldCell<Vect> {
  auto& m = eb.GetMesh();
  FieldCell<Vect> fcg(eb, Vect(0));
  for (auto c : eb.AllCells()) {
    if (eb.GetType(c) == Type::regular) {
      Vect sum(0);
      for (auto q : m.Nci(c)) {
        const IdxFace f = eb.GetFace(c, q);
        sum += m.GetOutwardSurface(c, q) * feu[f];
      }
      fcg[c] = sum / m.GetVolume(c);
    } else { // Type::cut
      std::vector<Vect> xx;
      std::vector<Scal> uu;
      xx.push_back(eb.GetFaceCenter(c));
      uu.push_back(feu[c]);
      for (auto q : eb.Nci(c)) {
        const IdxFace f = eb.GetFace(c, q);
        xx.push_back(eb.GetFaceCenter(f));
        uu.push_back(feu[f]);
      }
      auto p = ULinear<Scal>::FitLinear(xx, uu);
      fcg[c] = p.first;
    }
  }
  return fcg;
}

template <class M>
template <class T>
auto UEmbed<M>::AverageCutCells(const FieldCell<T>& fcu, const EB& eb)
    -> FieldCell<T> {
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

template <class M>
template <class T>
auto UEmbed<M>::RedistributeCutCells(const FieldCell<T>& fcu, const EB& eb)
    -> FieldCell<T> {
  auto& m = eb.GetMesh();
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
        fcr[c] += dun * (v / (eb.GetVolumeStencilSum(cn) - vn));
      }
    }
  }
  return fcr;
}

template <class M>
template <class T>
auto UEmbed<M>::Gradient(
    const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
    const EB& eb) -> FieldEmbed<T> {
  FieldEmbed<T> feg(eb, T(0));
  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    const Scal dn = eb.ClipGradDenom(
        (eb.GetNormal(f)).dot(eb.GetCellCenter(cp) - eb.GetCellCenter(cm)));
    feg[f] = (fcu[cp] - fcu[cm]) / dn;
  }
  for (auto c : eb.SuCFaces()) {
    if (bc == 0) {
      const Scal dn = eb.ClipGradDenom(eb.GetFaceOffset(c));
      feg[c] = (bcv - fcu[c]) / dn;
    } else if (bc == 1) {
      feg[c] = bcv;
    } else {
      throw std::runtime_error("Gradient: unknown bc=" + std::to_string(bc));
    }
  }
  GradientB(fcu, mfc, eb.GetMesh(), feg.GetFieldFace());
  return feg;
}

template <class M>
template <class T>
auto UEmbed<M>::InterpolateBilinearFaces(const FieldFace<T>& ffu, const EB& eb)
    -> FieldFace<T> {
  auto& m = eb.GetMesh();
  FieldFace<T> ffb(eb, T(0));
  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    const Scal h = eb.GetCellSize()[0];
    if (eb.GetType(cm) == Type::regular && eb.GetType(cp) == Type::regular) {
      ffb[f] = ffu[f];
    } else {
      const size_t qz = eb.GetNci(cm, f); // f == GetFace(c, qz)
      //                                             //
      //            ---------------------            //
      //            |         |         |            //
      //            |   c01   |   c11   |            //
      //            |         |         |            //
      //            |\--------|---------|            //
      //            | \       |         |            //
      //            |  \ c00  |   c10   |            //
      //            |   \     |         |            //
      //            -----\---------------            //
      //                                             //
      const size_t dz = qz / 2; // face direction
      const size_t dx = (dz + 1) % dim; // other directions
      const size_t dy = (dz + 2) % dim; // other directions
      const Vect n = eb.GetNormal(cm);
      const size_t qx = dx * 2 + (n[dx] < 0 ? 1 : 0);
      const size_t qy = dy * 2 + (n[dy] < 0 ? 1 : 0);
      const IdxCell c00 = cm;
      const IdxCell c10 = eb.GetCell(c00, qx);
      const IdxCell c01 = eb.GetCell(c00, qy);
      const IdxCell c11 = eb.GetCell(c10, qy);
      const IdxFace f00 = eb.GetFace(c00, qz);
      const IdxFace f10 = eb.GetFace(c10, qz);
      const IdxFace f01 = eb.GetFace(c01, qz);
      const IdxFace f11 = eb.GetFace(c11, qz);

      assert(f == f00);
      const Scal tx =
          1 - std::abs(eb.GetFaceCenter(f)[dx] - m.GetCenter(f)[dx]) / h;
      const Scal ty =
          1 - std::abs(eb.GetFaceCenter(f)[dy] - m.GetCenter(f)[dy]) / h;
      const T fy0 = ffu[f00] * tx + ffu[f10] * (1 - tx);
      const T fy1 = ffu[f01] * tx + ffu[f11] * (1 - tx);
      ffb[f] = fy0 * ty + fy1 * (1 - ty);
    }
  }
  return ffb;
}

template <class M>
template <class T>
auto UEmbed<M>::InterpolateBilinear(const FieldCell<T>& fcu, const EB& eb)
    -> FieldFace<T> {
  FieldFace<T> ffu(eb, T(0));
  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    ffu[f] = (fcu[cp] + fcu[cm]) * 0.5;
  }
  return InterpolateBilinearFaces(ffu, eb);
}

template <class M>
template <class T>
void UEmbed<M>::InterpolateEmbedFaces(
    const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu,
    FieldEmbed<T>& feu, const EB& eb) {
  auto& m = eb.GetMesh();
  if (bc == 0) {
    for (auto c : eb.SuCFaces()) {
      feu[c] = mcu.at(c);
    }
  } else if (bc == 1) {
    for (auto c : eb.SuCFaces()) {
      std::vector<Vect> xx;
      std::vector<T> uu;
      for (auto cn : eb.Stencil(c)) {
        xx.push_back(m.GetCenter(cn));
        uu.push_back(fcu[cn]);
      }
      auto p = ULinear<Scal>::FitLinear(xx, uu);
      const Scal h = eb.GetCellSize()[0];
      const Vect xc = eb.GetFaceCenter(c) - eb.GetNormal(c) * h;
      const T uc = ULinear<Scal>::EvalLinear(p, xc);
      feu[c] = uc + mcu.at(c) * h;
    }
  } else {
    throw std::runtime_error("not implemented");
  }
}

template <class M>
template <class T>
void UEmbed<M>::InterpolateUpwindEmbedFaces(
    const FieldCell<T>& fcu, size_t bc, const MapCell<Scal>& mcu,
    const FieldEmbed<Scal>& fev, FieldEmbed<T>& feu, const EB& eb) {
  auto& m = eb.GetMesh();
  if (bc == 0) {
    for (auto c : eb.SuCFaces()) {
      if (fev[c] > 0) {
        feu[c] = fcu[c];
      } else if (fev[c] <= 0) {
        feu[c] = mcu.at(c);
      }
    }
  } else if (bc == 1) {
    for (auto c : eb.SuCFaces()) {
      if (fev[c] > 0) {
        feu[c] = fcu[c];
      } else if (fev[c] <= 0) {
        std::vector<Vect> xx;
        std::vector<Scal> uu;
        for (auto cn : eb.Stencil(c)) {
          xx.push_back(m.GetCenter(cn));
          uu.push_back(fcu[cn]);
        }
        auto p = ULinear<Scal>::FitLinear(xx, uu);
        const Scal h = eb.GetCellSize()[0];
        const Vect xc = eb.GetFaceCenter(c) - eb.GetNormal(c) * h;
        const Scal uc = ULinear<Scal>::EvalLinear(p, xc);
        feu[c] = uc + mcu.at(c) * h;
      }
    }
  } else {
    throw std::runtime_error("not implemented");
  }
}

template <class M>
template <class T>
auto UEmbed<M>::InterpolateBilinear(
    const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu, const EB& eb)
    -> FieldEmbed<T> {
  FieldEmbed<T> feu(eb, T(0));
  feu.GetFieldFace() = InterpolateBilinear(fcu, eb);
  InterpolateEmbedFaces(fcu, bc, mcu, feu, eb);
  return feu;
}

template <class M>
template <class T>
auto UEmbed<M>::InterpolateBilinear(
    const FieldCell<T>& fcu, size_t bc, T bcv, const EB& eb) -> FieldEmbed<T> {
  MapCell<T> mcu;
  for (auto c : eb.SuCFaces()) {
    mcu[c] = bcv;
  }
  return InterpolateBilinear(fcu, bc, mcu, eb);
}

template <class M>
template <class T>
auto UEmbed<M>::InterpolateBilinear(
    const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
    const EB& eb) -> FieldEmbed<T> {
  MapCell<T> mcu;
  for (auto c : eb.SuCFaces()) {
    mcu[c] = bcv;
  }
  auto feu = InterpolateBilinear(fcu, bc, mcu, eb);
  InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
  return feu;
}

template <class M>
template <class T>
auto UEmbed<M>::GradientBilinear(const FieldCell<T>& fcu, const EB& eb)
    -> FieldFace<T> {
  FieldFace<T> ffg(eb, T(0));
  const Scal h = eb.GetCellSize()[0];
  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    ffg[f] = (fcu[cp] - fcu[cm]) / h;
  }
  return InterpolateBilinearFaces(ffg, eb);
}

template <class M>
template <class T>
auto UEmbed<M>::GradientBilinear(
    const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu, const EB& eb)
    -> FieldEmbed<T> {
  auto& m = eb.GetMesh();
  FieldEmbed<T> feg(eb, T(0));
  feg.GetFieldFace() = GradientBilinear(fcu, eb);
  if (bc == 0) {
    for (auto c : eb.SuCFaces()) {
      std::vector<Vect> xx;
      std::vector<T> uu;
      for (auto cn : eb.Stencil(c)) {
        xx.push_back(m.GetCenter(cn));
        uu.push_back(fcu[cn]);
      }
      auto p = ULinear<Scal>::FitLinear(xx, uu);
      const T ub = mcu.at(c);
      const Scal h = eb.GetCellSize()[0];
      const Vect xc = eb.GetFaceCenter(c) - eb.GetNormal(c) * h;
      const T uc = ULinear<Scal>::EvalLinear(p, xc);
      feg[c] = (ub - uc) / h;
    }
  } else if (bc == 1) {
    for (auto c : eb.SuCFaces()) {
      feg[c] = mcu.at(c);
    }
  } else {
    throw std::runtime_error("not implemented");
  }
  return feg;
}

template <class M>
template <class T>
auto UEmbed<M>::GradientBilinear(
    const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
    const EB& eb) -> FieldEmbed<T> {
  MapCell<T> mcu;
  for (auto c : eb.SuCFaces()) {
    mcu[c] = bcv;
  }
  auto feg = GradientBilinear(fcu, bc, mcu, eb);
  GradientB(fcu, mfc, eb.GetMesh(), feg.GetFieldFace());
  return feg;
}
