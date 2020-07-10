// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "approx_eb.h"

#include <inside.h>
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

// Returns level-set positive inside the body normalized to a new bounding box.
// center: center of the new bounding box
// extent: extent of the new bounding box
// rotation: rotation vector
// path: path to model
template <class M>
FieldNode<typename M::Scal> GetModelLevelSet(
    typename M::Vect center, typename M::Scal extent, typename M::Vect rotation,
    std::string path, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  int nt;
  int nv;
  int* tri;
  double* ver;
  struct Inside* inside_state;

  if (inside_mesh_read(path.c_str(), &nt, &tri, &nv, &ver) != 0) {
    throw std::runtime_error(FILELINE + ": fail to read mesh '" + path + "'");
  }

  // normalize ver to [0,1] and rotate
  {
    // rotates x around omega by angle |omega| degrees with origin at 0
    auto rotate = [](Vect x, Vect omega) {
      const Scal magn = omega.norm();
      if (magn == 0) {
        return x;
      }
      const Vect n = omega / magn;
      const Vect tu = Vect::GetUnit(omega.abs().argmin()).cross(n);
      const Vect tv = n.cross(tu);
      const Scal xu = x.dot(tu);
      const Scal xv = x.dot(tv);
      const Scal a = magn / 180 * M_PI;
      const Scal ru = xu * std::cos(a) - xv * std::sin(a);
      const Scal rv = xu * std::sin(a) + xv * std::cos(a);
      return tu * ru + tv * rv + n * x.dot(n);
    };
    // bounding box
    Vect box0(std::numeric_limits<Scal>::max());
    Vect box1(-std::numeric_limits<Scal>::max());
    for (int i = 0; i < nv * 3; ++i) {
      const int d = i % 3;
      box0[d] = std::min<Scal>(box0[d], ver[i]);
      box1[d] = std::max<Scal>(box1[d], ver[i]);
    }
    const Scal box_extent = (box1 - box0).abs().max();
    const Vect box_center = (box0 + box1) * 0.5;
    for (int i = 0; i < nv * 3; i += 3) {
      Vect x(ver[i], ver[i + 1], ver[i + 2]);
      x = center + rotate(x - box_center, rotation) / box_extent * extent;
      ver[i] = x[0];
      ver[i + 1] = x[1];
      ver[i + 2] = x[2];
    }
  }

  inside_ini(nt, tri, ver, &inside_state);
  auto distance = [&](Vect x) -> Scal {
    double p[3] = {x[0], x[1], x[2]};
    return inside_distance(inside_state, p);
  };
  auto inside = [&](Vect x) -> bool {
    double p[3] = {x[0], x[1], x[2]};
    return inside_inside(inside_state, p);
  };

  FieldNode<Scal> fnl(m, GetNan<Scal>());
  for (auto c : m.AllCells()) {
    const size_t num_nodes = m.GetNumNodes(c);
    size_t num_inside = 0;
    for (size_t i = 0; i < num_nodes; ++i) {
      const IdxNode n = m.GetNode(c, i);
      if (inside(m.GetNode(n))) {
        ++num_inside;
      }
    }
    // levelset > 0 inside body
    const Scal inf = m.GetCellSize()[0] * 10; // large value
    if (num_inside == 0) { // whole cell outside
      for (size_t i = 0; i < num_nodes; ++i) {
        const IdxNode n = m.GetNode(c, i);
        if (IsNan(fnl[n])) {
          fnl[n] = -inf;
        }
      }
    } else if (num_inside == num_nodes) { // whole cell inside
      for (size_t i = 0; i < num_nodes; ++i) {
        const IdxNode n = m.GetNode(c, i);
        if (IsNan(fnl[n])) {
          fnl[n] = inf;
        }
      }
    } else { // cell crosses surface
      for (size_t i = 0; i < num_nodes; ++i) {
        const IdxNode n = m.GetNode(c, i);
        fnl[n] = -distance(m.GetNode(n));
      }
    }
  }
  inside_fin(inside_state);
  inside_mesh_fin(tri, ver);
  return fnl;
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
  } else if (name == "model") {
    const std::string path = var.String["eb_model_path"];
    const Vect center(var.Vect["eb_model_center"]);
    const Scal extent = var.Double["eb_model_extent"];
    const Vect rotation(var.Vect["eb_model_rotation"]);
    fnl = GetModelLevelSet(center, extent, rotation, path, m);
  } else if (name == "list") {
    // TODO revise with bcast
    std::stringstream in;
    {
      std::stringstream path(var.String["eb_list_path"]);
      std::string filename;
      path >> filename;
      if (filename == "inline") {
        if (verb) {
          std::cout
              << "InitEmbed: Reading inline list of primitives from list_path"
              << std::endl;
        }
        in << path.rdbuf();
      } else {
        filename = path.str();
        std::ifstream fin(filename);
        if (verb) {
          std::cout << "InitEmbed: Open list of primitives '" << filename
                    << std::endl;
        }
        if (!fin.good()) {
          throw std::runtime_error(
              FILELINE + ": Can't open list of primitives");
        }
        in << fin.rdbuf();
      }
    }

    const size_t edim = var.Int["dim"];
    auto pp = UPrimList<Scal>::Parse(in, verb, edim);

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
    throw std::runtime_error(FILELINE + ": Unknown eb_init=" + name);
  }

  if (var.Int["eb_init_inverse"]) {
    for (auto n : m.AllNodes()) {
      fnl[n] = -fnl[n];
    }
  }
  fnl.SetHalo(2);
  return fnl;
}

// Evaluates quadratic interpolant on points -1, 0, 1.
// x: target point
// um,u,up: values of function for x=-1, 0, 1.
template <class Scal, class T>
T Quad(Scal x, T um, T u, T up) {
  return (um * (x - 1) + up * (x + 1)) * x * 0.5 - u * (x - 1) * (x + 1);
}


// Computes gradient at embeded face with Dirichlet conditions.
// rf: face center
// uf: given field value from boundary conditions
// nf: normal to face, towards excluded domain
// c: cell containign embedded face
// fcu: input field
template <class M, class T>
T GradDirichletQuad(
    typename M::Vect rf, T uf, typename M::Vect nf, IdxCell c,
    const FieldCell<T>& fcu, const M& m) {
  //            ---------------------            //
  //            |         |         |            //
  //            |         |  cw+dy  |            //
  //            |         |         |            //
  //            |\--------|---------|            //
  //            | \       |         |            //
  //            |  f  c   | cw=c+dx |            //
  //            |   \     |         |            //
  //            |----\----|---------|            //
  //            |     \   |         |            //
  //            |      \  |  cw-dy  |            //
  //            |       \ |         |            //
  //            |---------|---------|            //
  //                                             //
  auto h = m.GetCellSize()[0];
  auto dx = m.direction(nf.abs().argmax());
  auto dy = m.direction(1 - nf.abs().argmax()).orient(-nf);
  auto cw = m(c) + dx;
  auto dist = (cw.center[dx] - rf[dx]) / nf[dx];
  auto yi = rf[dy] + dist * nf[dy];
  auto ui = Quad((yi - cw.center[dy]) / h, fcu[cw - dy], fcu[cw], fcu[cw + dy]);
  return (uf - ui) / std::abs(dist);
}

template <class M, class T>
T GradDirichletLinear(
    typename M::Vect rf, T uf, typename M::Vect nf, IdxCell c,
    const FieldCell<T>& fcu, const Embed<M>& eb) {
  using Scal = typename M::Scal;
  auto linear = [](Scal x, T a0, T a1) { return a1 * x + a0 * (1 - x); };

  const auto& m = eb.GetMesh();
  const auto h = eb.GetCellSize()[0];

  const size_t dx = nf.abs().argmax();
  const size_t dy = (dx == 0 ? 1 : 0);
  const size_t sx = (nf[dx] > 0 ? 0 : 1); // step against normal
  const size_t sy = (nf[dy] > 0 ? 0 : 1);
  const IdxCell c0 = c;
  const IdxCell c1 = eb.GetCell(c0, 2 * dx + sx);
  const IdxCell c1p = eb.GetCell(c1, 2 * dy + sy);
  const Scal x1 = m.GetCenter(c1)[dx];
  const Scal y1 = m.GetCenter(c1)[dy];
  const Scal xf = rf[dx];
  const Scal yf = rf[dy];
  const Scal xf1 = x1;
  const Scal yf1 = yf + (xf1 - xf) / nf[dx] * abs(nf[dy]);
  const T u1 = linear(-(yf1 - y1) / h, fcu[c1], fcu[c1p]);
  return (uf - u1) / ((xf - x1) / nf[dx]);
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

  auto calc = [&](auto cf, IdxCell c, const BCond<T>& bc) {
    const Scal h = eb.GetCellSize()[0];
    const T& val = bc.val;
    switch (bc.type) {
      case BCondType::dirichlet: {
        return val;
      }
      case BCondType::neumann: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const T u = EvalLinearFit(x, c, fcu, eb);
        return u + val * h;
      }
      case BCondType::mixed:
      case BCondType::reflect: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const T u = EvalLinearFit(x, c, fcu, eb);
        return CombineMixed<Scal>()(val, u + val * h, eb.GetNormal(cf));
      }
      case BCondType::extrap: {
        return EvalLinearFit(eb.GetFaceCenter(cf), c, fcu, eb);
      }
    }
    return GetNan<T>();
  };
  mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto bc) { //
    feu[cf] = calc(cf, c, bc);
  });
  return feu;
}

template <class M>
template <class T>
auto UEmbed<M>::Gradient(
    const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const EB& eb)
    -> FieldEmbed<T> {
  FieldEmbed<T> feg(eb, T(0));
  const Scal h = eb.GetCellSize()[0];

  for (auto f : eb.SuFaces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    feg[f] = (fcu[cp] - fcu[cm]) / h;
  }
  feg.GetFieldFace() = InterpolateBilinearFaces(feg.GetFieldFace(), eb);

  auto calc = [&](auto cf, IdxCell c, const BCond<T>& bc) {
    const T& val = bc.val;
    switch (bc.type) {
      case BCondType::dirichlet: {
        //const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        //const T u = EvalLinearFit(x, c, fcu, eb);
        //return (val - u) / h;
        return GradDirichletQuad(
            eb.GetFaceCenter(cf), val, eb.GetNormal(cf), c, fcu, eb.GetMesh());
        // return GradDirichletLinear(
        //    eb.GetFaceCenter(cf), val, eb.GetNormal(cf), c, fcu, eb);
      }
      case BCondType::neumann: {
        return val;
      }
      case BCondType::mixed:
      case BCondType::reflect: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const T u = EvalLinearFit(x, c, fcu, eb);
        return CombineMixed<Scal>()((val - u) / h, val, eb.GetNormal(cf));
      }
      case BCondType::extrap: {
        // TODO replace with dot-product of gradient and normal
        auto p = FitLinear(c, fcu, eb);
        const Vect x1 = eb.GetFaceCenter(cf);
        const Vect x0 = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const T u1 = ULinear<typename M::Scal>::EvalLinear(p, x1);
        const T u0 = ULinear<typename M::Scal>::EvalLinear(p, x0);
        return (u1 - u0) / h;
      }
    }
    return GetNan<T>();
  };
  mebc.LoopBCond(
      eb, [&](auto cf, IdxCell c, auto bc) { feg[cf] = calc(cf, c, bc); });
  feg.LimitHalo(1);
  feg.LimitHalo(fcu.GetHalo() - 1);
  return feg;
}

template <class M>
auto UEmbed<M>::InterpolateUpwind(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    ConvSc scheme, const FieldCell<Vect>& fcg, const FieldEmbed<Scal>& fev,
    const EB& eb) -> FieldEmbed<Scal> {
  auto& m = eb.GetMesh();
  FieldEmbed<Scal> feu(eb, 0);
  if (scheme == ConvSc::superbee) {
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      const Scal du = fcu[cp] - fcu[cm];
      if (fev[f] > 0) {
        const auto rm = m.GetVectToCell(f, 0);
        feu[f] = fcu[cm] + 0.5 * Superbee(du, -4. * fcg[cm].dot(rm) - du);
      } else if (fev[f] < 0) {
        const auto rp = m.GetVectToCell(f, 1);
        feu[f] = fcu[cp] - 0.5 * Superbee(du, 4. * fcg[cp].dot(rp) - du);
      } else {
        feu[f] = (fcu[cm] + fcu[cp]) * 0.5;
      }
    }
  } else {
    const std::array<Scal, 3> a = GetCoeff<Scal>(scheme);
    // f = fmm*a[0] + fm*a[1] + fp*a[2]
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      if (fev[f] > 0) {
        const auto rm = m.GetVectToCell(f, 0);
        feu[f] = 4. * a[0] * fcg[cm].dot(rm) //
                 + a[1] * fcu[cm] + (a[2] + a[0]) * fcu[cp];
      } else if (fev[f] < 0) {
        const auto rp = m.GetVectToCell(f, 1);
        feu[f] = 4. * a[0] * fcg[cp].dot(rp) //
                 + a[1] * fcu[cp] + (a[2] + a[0]) * fcu[cm];
      } else {
        feu[f] = (fcu[cm] + fcu[cp]) * 0.5;
      }
    }
  }
  feu.GetFieldFace() = InterpolateBilinearFaces(feu.GetFieldFace(), eb);

  auto calc = [&](auto cf, IdxCell c, const BCond<Scal>& bc) {
    // TODO reuse code from Interpolate()
    const Scal h = eb.GetCellSize()[0];
    const Scal& val = bc.val;
    switch (bc.type) {
      case BCondType::dirichlet: {
        return val;
      }
      case BCondType::neumann: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const Scal u = EvalLinearFit(x, c, fcu, eb);
        return u + val * h;
      }
      case BCondType::mixed:
      case BCondType::reflect: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const Scal u = EvalLinearFit(x, c, fcu, eb);
        return CombineMixed<Scal>()(val, u + val * h, eb.GetNormal(cf));
      }
      case BCondType::extrap: {
        return EvalLinearFit(eb.GetFaceCenter(cf), c, fcu, eb);
      }
    }
    return GetNan<Scal>();
  };
  mebc.LoopBCond(eb, [&](auto cf, IdxCell c, const auto& bc) {
    feu[cf] = calc(cf, c, bc);
  });
  return feu;
}

template <class M>
auto UEmbed<M>::InterpolateBcg(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    const FieldEmbed<Scal>& fev, const FieldCell<Scal>& fc_src, const Scal dt,
    const EB& eb) -> FieldEmbed<Scal> {
  const Scal h = eb.GetCellSize()[0];

  FieldEmbed<bool> is_boundary(eb, false);
  mebc.LoopPairs([&](auto p) { //
    is_boundary[p.first] = true;
  });
  for (auto f : eb.SuFaces()) {
    if (eb.GetType(f) != EB::Type::regular) {
      is_boundary[f] = true;
    }
  }

  const FieldEmbed<Scal> feg = Gradient(fcu, mebc, eb);
  FieldEmbed<Scal> feu(eb, 0);
  for (auto f : eb.Faces()) {
    const Scal sgn = (fev[f] > 0 ? 1 : -1);
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    const IdxCell c = (sgn > 0 ? cm : cp); // upwind cell
    const size_t nci = eb.GetNci(c, f);
    auto Fm = [&](size_t d) { //
      return eb.GetFace(c, ((nci / 2 + d) % 3) * 2);
    };
    auto Fp = [&](size_t d) { //
      return eb.GetFace(c, ((nci / 2 + d) % 3) * 2 + 1);
    };

    // temporal derivative, initial from acceleration
    Scal ut = (fc_src[cm] + fc_src[cp]) * 0.5;
    // convective terms
    // normal direction
    const IdxFace fm = Fm(0);
    const IdxFace fp = Fp(0);
    const Scal ux = (is_boundary[fm] || is_boundary[fp])
                        ? feg[f]
                        : (feg[fm] + feg[fp]) * 0.5; // normal derivative
    const Scal w = fev[f] / eb.GetArea(f); // velocity
    ut -= ux * w;
    // tangential directions
    for (auto d : {1, 2}) {
      const IdxFace fm = Fm(d);
      const IdxFace fp = Fp(d);
      if (is_boundary[fm] || is_boundary[fp]) {
        // exclude faces that depend on boundary conditions
        // as they cause instability in case of non-zero Dirichlet conditions
        continue;
      }
      const Scal w = (fev[fm] + fev[fp]) / (eb.GetArea(fm) + eb.GetArea(fp));
      ut -= feg[w > 0 ? fm : fp] * w;
    }

    feu[f] = fcu[c] + ux * (sgn * 0.5 * h) + ut * (0.5 * dt);
  }
  auto calc = [&](auto cf, IdxCell c, const BCond<Scal>& bc) {
    const Scal& val = bc.val;
    switch (bc.type) {
      case BCondType::dirichlet: {
        return val;
      }
      case BCondType::neumann: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const Scal u = EvalLinearFit(x, c, fcu, eb);
        return u + val * h;
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    return GetNan<Scal>();
  };
  mebc.LoopBCond(eb, [&](auto cf, IdxCell c, const auto& bc) {
    feu[cf] = calc(cf, c, bc);
  });
  feu.LimitHalo(0);
  feu.LimitHalo(fcu.GetHalo() - 1);
  feu.LimitHalo(fc_src.GetHalo() - 1);
  feu.LimitHalo(fev.GetHalo() - 1);
  return feu;
}

template <class M>
auto UEmbed<M>::InterpolateBcg(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    const FieldFace<Scal>& ffv, const FieldCell<Scal>& fc_src, const Scal dt,
    const M& m) -> FieldFace<Scal> {
  const Scal h = m.GetCellSize()[0];

  FieldEmbed<bool> is_boundary(m, false);
  mebc.LoopPairs([&](auto p) { //
    is_boundary[p.first] = true;
  });

  const FieldFace<Scal> ffg = Gradient(fcu, mebc, m);
  FieldFace<Scal> ffu(m, 0);
  for (auto f : m.Faces()) {
    const Scal sgn = (ffv[f] > 0 ? 1 : -1);
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    const IdxCell c = (sgn > 0 ? cm : cp); // upwind cell
    const size_t nci = m.GetNci(c, f);
    auto Fm = [&](size_t d) { //
      return m.GetFace(c, ((nci / 2 + d) % 3) * 2);
    };
    auto Fp = [&](size_t d) { //
      return m.GetFace(c, ((nci / 2 + d) % 3) * 2 + 1);
    };

    // temporal derivative, initial from acceleration
    Scal ut = (fc_src[cm] + fc_src[cp]) * 0.5;
    // convective terms
    // normal direction
    const IdxFace fm = Fm(0);
    const IdxFace fp = Fp(0);
    const Scal ux = (ffg[fm] + ffg[fp]) * 0.5; // normal derivative
    const Scal w = ffv[f] / m.GetArea(f); // velocity
    ut -= ux * w;
    // tangential directions
    for (auto d : {1, 2}) {
      const IdxFace fm = Fm(d);
      const IdxFace fp = Fp(d);
      if (is_boundary[fm] || is_boundary[fp]) {
        // exclude faces that depend on boundary conditions
        // as they cause instability in case of non-zero Dirichlet conditions
        continue;
      }
      const Scal w = (ffv[fm] + ffv[fp]) / (m.GetArea(fm) + m.GetArea(fp));
      ut -= ffg[w > 0 ? fm : fp] * w;
    }

    ffu[f] = fcu[c] + ux * (sgn * 0.5 * h) + ut * (0.5 * dt);
  }
  auto calc = [&](IdxFace f, IdxCell c, const BCond<Scal>& bc) {
    const Scal& val = bc.val;
    const auto nci = bc.nci;
    switch (bc.type) {
      case BCondType::dirichlet: {
        return val;
      }
      case BCondType::neumann: {
        const Scal q = (nci == 0 ? 1. : -1.);
        const Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * q;
        return fcu[c] + bc.val * a;
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    return GetNan<Scal>();
  };
  // faces with boundary conditions
  for (auto& p : mebc.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    ffu[f] = calc(f, m.GetCell(f, bc.nci), bc);
  }
  ffu.LimitHalo(0);
  ffu.LimitHalo(fcu.GetHalo() - 1);
  ffu.LimitHalo(fc_src.GetHalo() - 1);
  ffu.LimitHalo(ffv.GetHalo() - 1);
  return ffu;
}

template <class M>
template <class T>
auto UEmbed<M>::Interpolate(
    const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const M& m)
    -> FieldFace<T> {
  FieldFace<T> ffu(m, T(0));

  for (auto f : m.SuFaces()) {
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    ffu[f] = (fcu[cp] + fcu[cm]) * 0.5;
  }

  auto calc = [&](IdxFace f, IdxCell c, const BCond<T>& bc) {
    const T& val = bc.val;
    const auto nci = bc.nci;
    switch (bc.type) {
      case BCondType::dirichlet: {
        return val;
      }
      case BCondType::neumann: {
        const Scal q = (nci == 0 ? 1. : -1.);
        const Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * q;
        return fcu[c] + val * a;
      }
      case BCondType::mixed:
      case BCondType::reflect: {
        const Scal q = (nci == 0 ? 1. : -1.);
        const Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * q;
        return CombineMixed<Scal>()(val, fcu[c] + val * a, m.GetNormal(f));
      }
      case BCondType::extrap: {
        const IdxCell c = m.GetCell(f, nci);
        const size_t q = m.GetNci(c, f);
        const size_t qo = m.GetOpposite(q);
        const IdxFace fo = m.GetFace(c, qo);
        const Vect n = m.GetNormal(f);
        // cell
        const T& v0 = fcu[c];
        const Scal x0 = 0.;
        // opposite face
        const T& v1 = ffu[fo];
        const Scal x1 = n.dot(m.GetCenter(fo) - m.GetCenter(c));
        // target
        const Scal xt = n.dot(m.GetCenter(f) - m.GetCenter(c));
        return UExtrap(xt, x0, v0, x1, v1);
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    return GetNan<T>();
  };
  // faces with boundary conditions
  for (auto& p : mebc.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    ffu[f] = calc(f, m.GetCell(f, bc.nci), bc);
  }
  return ffu;
}

template <class M>
template <class T>
auto UEmbed<M>::Gradient(
    const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const M& m)
    -> FieldFace<T> {
  FieldFace<T> ffu(m, T(0));
  const Scal h = m.GetCellSize()[0];

  for (auto f : m.SuFaces()) {
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    ffu[f] = (fcu[cp] - fcu[cm]) / h;
  }

  auto calc = [&](IdxFace f, IdxCell c, const BCond<T>& bc) {
    const T& val = bc.val;
    const auto nci = bc.nci;
    switch (bc.type) {
      case BCondType::dirichlet: {
        const Scal sgn = (nci == 0 ? 1 : -1);
        const T uf = val;
        const T u1 = fcu[c];
        const size_t cnci = m.GetNci(c, f); // f=m.GetFace(c, cnci)
        const T u2 = fcu[m.GetCell(c, m.GetOpposite(cnci))];
        //    |      |      |
        //  uf|  u1  |  u2  |
        //    |      |      |
        return (uf * 8 - u1 * 9 + u2) / (3 * h) * sgn;
      }
      case BCondType::neumann: {
        return val;
      }
      case BCondType::extrap: {
        const Scal sgn = (nci == 0 ? 1 : -1);
        const T u1 = fcu[c];
        const size_t cnci = m.GetNci(c, f); // f=m.GetFace(c, cnci)
        const T u2 = fcu[m.GetCell(c, m.GetOpposite(cnci))];
        return (u1 - u2) / h * sgn;
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    return GetNan<T>();
  };
  // faces with boundary conditions
  for (auto& p : mebc.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    ffu[f] = calc(f, m.GetCell(f, bc.nci), bc);
  }
  ffu.LimitHalo(1);
  ffu.LimitHalo(fcu.GetHalo() - 2);
  return ffu;
}

template <class M>
auto UEmbed<M>::InterpolateUpwind(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    ConvSc scheme, const FieldCell<Vect>& fcg, const FieldFace<Scal>& ffv,
    const M& m) -> FieldFace<Scal> {
  FieldFace<Scal> ffu(m, 0);

  if (scheme == ConvSc::superbee) {
    for (auto f : m.SuFaces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const Scal du = fcu[cp] - fcu[cm];
      if (ffv[f] > 0) {
        const auto rm = m.GetVectToCell(f, 0);
        ffu[f] = fcu[cm] + 0.5 * Superbee(du, -4. * fcg[cm].dot(rm) - du);
      } else if (ffv[f] < 0) {
        const auto rp = m.GetVectToCell(f, 1);
        ffu[f] = fcu[cp] - 0.5 * Superbee(du, 4. * fcg[cp].dot(rp) - du);
      } else {
        ffu[f] = (fcu[cm] + fcu[cp]) * 0.5;
      }
    }
  } else {
    // f = fmm*a[0] + fm*a[1] + fp*a[2]
    std::array<Scal, 3> a = GetCoeff<Scal>(scheme);
    ffu.Reinit(m);
    for (auto f : m.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      if (ffv[f] > 0) {
        const auto rm = m.GetVectToCell(f, 0);
        ffu[f] = 4. * a[0] * fcg[cm].dot(rm) //
                 + a[1] * fcu[cm] + (a[2] + a[0]) * fcu[cp];
      } else if (ffv[f] < 0) {
        const auto rp = m.GetVectToCell(f, 1);
        ffu[f] = 4. * a[0] * fcg[cp].dot(rp) //
                 + a[1] * fcu[cp] + (a[2] + a[0]) * fcu[cm];
      } else {
        ffu[f] = (fcu[cm] + fcu[cp]) * 0.5;
      }
    }
  }

  auto calc = [&](IdxFace f, IdxCell c, const BCond<Scal>& bc) {
    const Scal& val = bc.val;
    const auto nci = bc.nci;
    switch (bc.type) {
      case BCondType::dirichlet: {
        return val;
      }
      case BCondType::neumann: {
        const Scal q = (nci == 0 ? 1. : -1.);
        const Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * q;
        return fcu[c] + bc.val * a;
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    return GetNan<Scal>();
  };
  // faces with boundary conditions
  for (auto& p : mebc.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    ffu[f] = calc(f, m.GetCell(f, bc.nci), bc);
  }
  return ffu;
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
auto UEmbed<M>::GradientGauss(const FieldEmbed<Scal>& feu, const EB& eb)
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
  for (auto c : eb.SuCells()) {
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
  fcu.CheckHalo(1);
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
auto UEmbed<M>::RedistributeCutCells(const FieldCell<T>& fcu, const M&)
    -> FieldCell<T> {
  return fcu;
}

// fcs: sum of fluxes
template <class M>
auto UEmbed<M>::RedistributeCutCellsAdvection(
    const FieldCell<Scal>& fcs, const FieldFace<Scal>&, Scal, Scal, const M&)
    -> FieldCell<Scal> {
  return fcs;
}

// fcs: sum of fluxes
template <class M>
auto UEmbed<M>::RedistributeCutCellsAdvection(
    const FieldCell<Scal>& fcs, const FieldEmbed<Scal>& ffv, Scal cfl, Scal dt,
    const Embed<M>& eb) -> FieldCell<Scal> {
  fcs.CheckHalo(1);
  // stability criterion:
  // dt < cfl * h / velocity
  // dt < cfl * V / volumeflux
  // dt * volumeflux < cfl * V
  // dtmax = cfl * V / volumeflux
  auto excess = [&](IdxCell c) {
    Scal dtmax = std::numeric_limits<Scal>::max();
    eb.LoopNci(c, [&](auto q) {
      const auto cf = eb.GetFace(c, q);
      const Scal flux = std::abs(ffv[cf]);
      if (flux != 0) {
        dtmax = std::min(dtmax, cfl * eb.GetVolume(c) / flux);
      }
    });
    return fcs[c] * std::max(0., (dt - dtmax) / dt);
  };
  auto fcr = fcs;
  for (auto c : eb.Cells()) {
    // full step:
    //   s * dt
    // allowed step:
    //   s * dtmax
    // excess step:
    //   s * (dt - dtmax)
    // excess sum of fluxes to be advected with dt:
    //   s * (dt - dtmax) / dt
    // excess sum of fluxes
    const Scal ds = excess(c);
    // subtract from current cell
    fcr[c] -= ds;
    // add from neighbor cells proportional to their volume
    for (auto cn : eb.Stencil(c)) {
      if (c != cn) {
        const Scal v = eb.GetVolume(c);
        const Scal vn = eb.GetVolume(cn);
        // excess quantity in cell cn
        const Scal dsn = excess(cn);
        fcr[c] += dsn * (v / (eb.GetVolumeStencilSum(cn) - vn));
      }
    }
  }
  return fcr;
}

template <class M>
void UEmbed<M>::RedistributeConstTerms(
    FieldCell<typename M::Expr>& fce, const Embed<M>& eb, M& m) {
  auto sem = m.GetSem("redistr");
  struct {
    FieldCell<Scal> fcu;
  } * ctx(sem);
  auto& fcu = ctx->fcu;
  auto get = [&fce](IdxCell c) -> Scal& { //
    return fce[c].back();
  };
  if (sem("comm")) {
    fcu.Reinit(eb, 0);
    for (auto c : eb.Cells()) {
      fcu[c] = get(c);
    }
    m.Comm(&fcu);
  }
  if (sem()) {
    // FIXME: empty stage to finish communication in halo cells
    // fcu gets different buffer in the next stage
  }
  if (sem("local")) {
    fcu = RedistributeCutCells(fcu, eb);
    for (auto c : eb.Cells()) {
      get(c) = fcu[c];
    }
  }
  if (sem()) {
    // FIXME: empty stage to finish communication and keep ctx
  }
}
template <class M>
void UEmbed<M>::RedistributeConstTerms(
    FieldCell<typename M::Expr>&, const M&, M&) {
  return;
}

template <class M>
template <class T>
auto UEmbed<M>::InterpolateBilinearFaces(const FieldFace<T>& ffu, const M&)
    -> FieldFace<T> {
  return ffu;
}

// Evaluates bilinear interpolant on points (0,0), (1,0), (0,1) and (1,1).
// x,y: target point
// u,ux,uy,uyx:  values of function for (x,y) = (0,0), (1,0), (0,1), (1,1)
template <class T, class Scal>
T Bilinear(Scal x, Scal y, T u, T ux, T uy, T uyx) {
  //                      //
  //   y                  //
  //   |                  //
  //   |*uy    *uyx       //
  //   |                  //
  //   |                  //
  //   |*u     *ux        //
  //   |-------------x    //
  //                      //
  const auto v = u * (1 - x) + ux * x;
  const auto vy = uy * (1 - x) + uyx * x;
  return v * (1 - y) + vy * y;
}

template <class M>
template <class T>
auto UEmbed<M>::InterpolateBilinearFaces(const FieldFace<T>& ffu, const EB& eb)
    -> FieldFace<T> {
  FieldFace<T> ffb(eb, T(0));
  for (auto f : eb.SuFacesM()) {
    const Scal h = eb.GetCellSize()[0];
    if (eb.GetType(f) == EB::Type::regular) {
      ffb[f] = ffu[f];
    } else {
      //   ---------------------  //
      //   |         |         |  //
      //   |    fy   |   fyx   |  //
      //   |         |         |  //
      //   |\--------|---------|  //
      //   | \       |         |  //
      //   |  r  f   |   fx    |  //
      //   |   \     |         |  //
      //   -----\---------------  //
      auto n = eb.GetNormal(f.cm);
      auto dz = f.direction();
      auto dx = (dz >> 1).orient(-n);
      auto dy = (dz >> 2).orient(-n);
      auto r = (eb.GetFaceCenter(f) - f.center).abs();
      ffb[f] = Bilinear(
          r[dx] / h, r[dy] / h, //
          ffu[f], ffu[f + dx], ffu[f + dy], ffu[f + dy + dx]);
    }
  }
  return ffb;
}

template <class M>
auto UEmbed<M>::InterpolateUpwindImplicit(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    ConvSc scheme, Scal deferred, const FieldCell<Vect>& fcg,
    const FieldEmbed<Scal>& fev, const EB& eb) -> FieldEmbed<ExprFace> {
  (void)deferred;
  FieldEmbed<ExprFace> fee(eb, ExprFace(0));

  // explicit interpolation
  FieldEmbed<Scal> feu = InterpolateUpwind(fcu, mebc, scheme, fcg, fev, eb);

  // implicit interpolation with deferred correction
  for (auto f : eb.Faces()) {
    ExprFace e(0);
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    if (fev[f] > 0) {
      e[0] = 1;
      e[2] = feu[f];
    } else if (fev[f] < 0) {
      e[1] = 1;
      e[2] = feu[f];
    } else {
      e[0] = 0.5;
      e[1] = 0.5;
      e[2] = feu[f];
    }
    e[2] -= fcu[cm] * e[0] + fcu[cp] * e[1];
    fee[f] = e;
  }

  auto calc = [&](auto cf, IdxCell c, const BCond<Scal>& bc) {
    ExprFace e(0);
    const auto nci = bc.nci;
    switch (bc.type) {
      case BCondType::dirichlet: {
        e[2] = feu[cf];
        break;
      }
      case BCondType::neumann: {
        e[nci] = 1;
        e[2] = feu[cf];
        break;
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    e[2] -= fcu[c] * e[nci];
    return e;
  };
  mebc.LoopBCond(eb, [&](auto cf, IdxCell c, const auto& bc) {
    fee[cf] = calc(cf, c, bc);
  });
  return fee;
}

template <class M>
auto UEmbed<M>::InterpolateUpwindImplicit(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    ConvSc scheme, Scal deferred, const FieldCell<Vect>& fcg,
    const FieldFace<Scal>& ffv, const M& m) -> FieldFace<ExprFace> {
  FieldFace<ExprFace> ffe(m, ExprFace(0));

  // f = fmm*a[0] + fm*a[1] + fp*a[2]
  std::array<Scal, 3> a = GetCoeff<Scal>(scheme);
  const Scal df = deferred;
  const Scal dfm = 1 - df;
  for (auto f : m.Faces()) {
    ExprFace e(0);
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    if (ffv[f] > 0) {
      e[0] = a[1] * dfm + df;
      e[1] = (a[2] + a[0]) * dfm;
      e[2] = 4. * a[0] * fcg[cm].dot(m.GetVectToCell(f, 0)) +
             (a[1] - 1) * df * fcu[cm] + (a[2] + a[0]) * df * fcu[cp];
    } else if (ffv[f] < 0) {
      e[0] = (a[2] + a[0]) * dfm;
      e[1] = a[1] * dfm + df;
      e[2] = 4. * a[0] * fcg[cp].dot(m.GetVectToCell(f, 1)) +
             (a[1] - 1) * df * fcu[cp] + (a[2] + a[0]) * df * fcu[cm];
    } else {
      e[0] = 0.5;
      e[1] = 0.5;
    }
    ffe[f] = e;
  }

  auto calc = [&](IdxFace f, IdxCell c, const BCond<Scal>& bc) {
    ExprFace e(0);
    const auto nci = bc.nci;
    switch (bc.type) {
      case BCondType::dirichlet: {
        e[2] = bc.val;
        break;
      }
      case BCondType::neumann: {
        const Scal g = (nci == 0 ? 1 : -1);
        const Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * g;
        e[nci] = 1;
        e[2] = a * bc.val;
        break;
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    return e;
  };
  // faces with boundary conditions
  for (auto& p : mebc.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    ffe[f] = calc(f, m.GetCell(f, bc.nci), bc);
  }
  return ffe;
}

template <class M>
auto UEmbed<M>::GradientImplicit(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, const EB& eb)
    -> FieldEmbed<ExprFace> {
  FieldEmbed<ExprFace> fee(eb, ExprFace(0));
  const Scal h = eb.GetCellSize()[0];

  // explicit gradient
  const FieldEmbed<Scal> feg = Gradient(fcu, mebc, eb);

  // implicit gradient with deferred correction
  for (auto f : eb.Faces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    const Scal a = 1 / h;
    ExprFace e(0);
    e[0] = -a;
    e[1] = a;
    e[2] = feg[f];
    e[2] -= fcu[cm] * e[0] + fcu[cp] * e[1];
    fee[f] = e;
  }

  auto calc = [&](auto cf, IdxCell c, const BCond<Scal>& bc) {
    ExprFace e(0);
    const auto nci = bc.nci;
    switch (bc.type) {
      case BCondType::dirichlet: {
        const Scal a = 2 * (nci == 0 ? 1 : -1) / h;
        e[nci] = -a;
        e[2] = feg[cf];
        break;
      }
      case BCondType::neumann: {
        e[2] = feg[cf];
        break;
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    e[2] -= fcu[c] * e[nci]; // subtract explicit part
    return e;
  };
  mebc.LoopBCond(eb, [&](auto cf, IdxCell c, const auto& bc) {
    fee[cf] = calc(cf, c, bc);
  });
  return fee;
}

template <class M>
auto UEmbed<M>::GradientImplicit(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, const M& m)
    -> FieldFace<ExprFace> {
  (void)fcu;
  FieldFace<ExprFace> ffe(m, ExprFace(0));
  const Scal h = m.GetCellSize()[0];

  for (auto f : m.Faces()) {
    ExprFace e(0);
    const Scal a = 1 / h;
    e[0] = -a;
    e[1] = a;
    ffe[f] = e;
  }

  auto calc = [&](IdxFace f, IdxCell c, const BCond<Scal>& bc) {
    ExprFace e(0);
    const auto nci = bc.nci;
    switch (bc.type) {
      case BCondType::dirichlet: {
        const Scal sgn = (nci == 0 ? 1 : -1);
        const Scal a = 2 * sgn / h;
        e[nci] = -a;
        e[2] = bc.val * a;

        const Scal uf = bc.val;
        const size_t cnci = m.GetNci(c, f); // f=m.GetFace(c1, c1nci)
        const IdxCell c2 = m.GetCell(c, m.GetOpposite(cnci));
        const Scal u1 = fcu[c];
        const Scal u2 = fcu[c2];
        const Scal high = (uf * 8 - u1 * 9 + u2) / (3 * h) * sgn;
        const Scal low = (uf - u1) * a;
        e[2] += high - low;
        break;
      }
      case BCondType::neumann: {
        e[2] = bc.val;
        break;
      }
      default:
        throw std::runtime_error(FILELINE + ": unknown");
    }
    return e;
  };
  // faces with boundary conditions
  for (auto& p : mebc.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    ffe[f] = calc(f, m.GetCell(f, bc.nci), bc);
  }
  return ffe;
}

template <class M>
auto UEmbed<M>::Gradient(const FieldFace<Scal>& ffu, const M& m)
    -> FieldCell<Vect> {
  FieldCell<Vect> fcg(m, Vect(0));
  for (auto c : m.SuCells()) {
    Vect sum(0);
    for (auto q : m.Nci(c)) {
      const IdxFace f = m.GetFace(c, q);
      sum += m.GetOutwardSurface(c, q) * ffu[f];
    }
    fcg[c] = sum / m.GetVolume(c);
  }
  return fcg;
}

// FIXME: rename to AverageVect
template <class M>
auto UEmbed<M>::AverageGradient(const FieldEmbed<Scal>& feg, const EB& eb)
    -> FieldCell<Vect> {
  FieldCell<Vect> fcg(eb, Vect(0));
  for (auto c : eb.AllCells()) {
    Vect sum(0);
    Vect sumw(0);
    for (auto q : eb.Nci(c)) {
      const auto f = eb.GetFace(c, q);
      const Vect w = eb.GetSurface(f);
      sum += w * feg[f];
      sumw += w;
    }
    fcg[c] = sum / sumw;
  }
  fcg.LimitHalo(2);
  fcg.LimitHalo(feg.GetHalo() - 1);
  return fcg;
}

template <class M>
auto UEmbed<M>::AverageGradient(const FieldFace<Scal>& ffg, const M& m)
    -> FieldCell<Vect> {
  FieldCell<Vect> fcg(m, Vect(0));
  for (auto c : m.AllCells()) {
    Vect sum(0);
    for (auto q : m.Nci(c)) {
      const auto f = m.GetFace(c, q);
      sum += m.GetNormal(f) * (ffg[f] * 0.5);
    }
    fcg[c] = sum;
  }
  fcg.LimitHalo(2);
  fcg.LimitHalo(ffg.GetHalo() - 1);
  return fcg;
}
