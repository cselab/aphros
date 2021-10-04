// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "approx_eb.h"

#include <iostream>

#include "approx.h"
#include "func/primlist.h"
#include "inside/inside.h"
#include "util/format.h"

template <class Vect_>
auto ULinearFit<Vect_>::FitLinear(
    const std::vector<Vect>& xx, const std::vector<Scal>& uu)
    -> std::pair<Vect, Scal> {
  assert(xx.size() == uu.size());
  // sum 0.5 * [ (g.dot(x[k]) + u0 - u[k]) ** 2 ] -> min
  using Int = size_t;
  static constexpr Int N = dim + 1;
  std::array<Scal, N * N> a;
  std::array<Scal, N> b;
  auto aa = [&a](Int i, Int j) -> Scal& { //
    return a[i * N + j];
  };
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
  return {Vect(&v[0]), v[dim]};
}

template <class Vect_>
auto ULinearFit<Vect_>::FitLinear(
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
    for (size_t i = 0; i < dim; ++i) {
      p.first[i][d] = pd.first[i];
    }
  }
  return p;
}

// Returns level-set positive inside the body normalized to a new bounding box.
// center: center of the new bounding box
// extent: extent of the new bounding box
// rotation: rotation vector
// path: path to model
template <class M>
void InitLevelSetFromModel(
    FieldNode<typename M::Scal>& fnl, typename M::Vect center,
    typename M::Scal extent, typename M::Vect rotation, std::string path,
    M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Vect3 = generic::Vect<Scal, 3>;
  auto sem = m.GetSem(__func__);
  struct Context {
    int nt;
    int nv;
    int* tri;
    double* ver;
    Inside* inside_state;
    std::vector<Context*> ctx_all;
  };
  Context* ctx(sem);
  auto& t = *ctx;
  if (sem("gather")) {
    t.ctx_all.push_back(ctx);
    m.GatherToLead(&t.ctx_all);
  }
  if (sem("ini")) {
    if (m.IsLead()) {
      if (inside_mesh_read(path.c_str(), &t.nt, &t.tri, &t.nv, &t.ver) != 0) {
        fassert(false, "failed to read mesh '" + path + "'");
      }

      // rescale, translate and rotate the vertices
      {
        // rotates x around axes 0,1,2
        // by angles omega[0], omega[1], omega[2] degrees with origin at 0
        auto rotate = [&m](Vect x, Vect omega) {
          for (auto d : m.dirs) {
            const Scal a = omega[d] / 180 * M_PI;
            const size_t du = (d + 1) % m.dim;
            const size_t dv = (d + 2) % m.dim;
            const Scal u = x[du] * std::cos(a) - x[dv] * std::sin(a);
            const Scal v = x[du] * std::sin(a) + x[dv] * std::cos(a);
            x[du] = u;
            x[dv] = v;
          }
          return x;
        };
        // bounding box
        Vect box0(std::numeric_limits<Scal>::max());
        Vect box1(-std::numeric_limits<Scal>::max());
        for (int i = 0; i < t.nv * 3; ++i) {
          const int d = i % 3;
          box0[d] = std::min<Scal>(box0[d], t.ver[i]);
          box1[d] = std::max<Scal>(box1[d], t.ver[i]);
        }
        const Scal box_extent = (box1 - box0).abs().max();
        const Vect box_center = (box0 + box1) * 0.5;
        for (int i = 0; i < t.nv * 3; i += 3) {
          Vect3 x(t.ver[i], t.ver[i + 1], t.ver[i + 2]);
          x = Vect3(
              center +
              rotate(Vect(x) - box_center, rotation) / box_extent * extent);
          t.ver[i] = x[0];
          t.ver[i + 1] = x[1];
          t.ver[i + 2] = x[2];
        }
      }

      inside_ini(t.nt, t.tri, t.ver, &t.inside_state);

      if (m.IsRoot()) {
        InsideInfo info;
        inside_info(t.inside_state, &info);
        fassert(
            info.size > m.GetCellSize()[0],
            util::Format(
                "{}: Required info.size={} > h={} for model '{}'. Use finer "
                "mesh or coarser model",
                __func__, info.size, m.GetCellSize()[0], path));
      }

      for (auto* other : t.ctx_all) {
        other->nt = t.nt;
        other->nv = t.nv;
        other->tri = t.tri;
        other->ver = t.ver;
        other->inside_state = t.inside_state;
      }
    }
  }
  if (sem("local")) {
    // only valid if `x` is close to the surface
    auto distance = [&](Vect3 x) -> Scal {
      double p[3] = {x[0], x[1], x[2]};
      return inside_distance(t.inside_state, p);
    };
    auto inside = [&](Vect3 x) -> bool {
      double p[3] = {x[0], x[1], x[2]};
      return inside_inside(t.inside_state, p);
    };

    const Scal inf = m.GetCellSize()[0] * 10; // large value
    fnl.Reinit(m, GetNan<Scal>());

    FieldNode<bool> fn_inside(m);
    for (auto n : m.AllNodes()) {
      fn_inside[n] = inside(Vect3(m.GetNode(n)));
    }

    for (auto c : m.AllCells()) {
      const size_t num_nodes = m.GetNumNodes(c);
      size_t num_inside = 0;
      for (size_t i = 0; i < num_nodes; ++i) {
        const IdxNode n = m.GetNode(c, i);
        if (fn_inside[n]) {
          ++num_inside;
        }
      }
      // levelset > 0 inside body
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
          fnl[n] = 0; // mark to call distance() later
        }
      }
    }
    for (auto n : m.AllNodes()) {
      if (fnl[n] == 0) {
        fnl[n] = -distance(Vect3(m.GetNode(n)));
      }
    }
  }
  if (sem("fin")) {
    if (m.IsLead()) {
      inside_fin(t.inside_state);
      inside_mesh_fin(t.tri, t.ver);
    }
  }
}

template <class M>
void UEmbed<M>::InitLevelSet(
    FieldNode<Scal>& fnl, M& m, const Vars& var, bool verb) {
  auto sem = m.GetSem(__func__);
  if (sem("init")) {
    fnl.Reinit(m);
  }
  const auto name = var.String["eb_init"];
  if (name == "none") {
    if (sem()) {
      fnl.Reinit(m, 1);
    }
  } else if (name == "box") {
    if (sem()) {
      const Vect xc(var.Vect["eb_box_c"]);
      const Vect r(var.Vect["eb_box_r"]);
      const Scal angle = M_PI * var.Double["eb_box_angle"];
      auto rotate = [angle](Vect x) {
        const Scal sin = std::sin(angle);
        const Scal cos = std::cos(angle);
        Vect res(x);
        res[0] = x[0] * cos - x[1] * sin;
        res[1] = x[0] * sin + x[1] * cos;
        return res;
      };
      for (auto n : m.AllNodes()) {
        fnl[n] = (1 - (rotate(m.GetNode(n) - xc) / r).norminf()) *
                 (r / m.GetCellSize()).min();
      }
    }
  } else if (name == "sphere") {
    if (sem()) {
      const Vect xc(var.Vect["eb_sphere_c"]);
      const Vect r(var.Vect["eb_sphere_r"]);
      const Scal angle = M_PI * var.Double["eb_sphere_angle"];
      auto rotate = [angle](Vect x) {
        const Scal sin = std::sin(angle);
        const Scal cos = std::cos(angle);
        Vect res(x);
        res[0] = x[0] * cos - x[1] * sin;
        res[1] = x[0] * sin + x[1] * cos;
        return res;
      };
      for (auto n : m.AllNodes()) {
        fnl[n] = (rotate(m.GetNode(n) - xc) / r).norm() - 1;
      }
    }
  } else if (name == "model") {
    if (sem.Nested()) {
      const std::string path = var.String["eb_model_path"];
      const Vect center(var.Vect["eb_model_center"]);
      const Scal extent = var.Double["eb_model_extent"];
      const Vect rotation(var.Vect["eb_model_rotation"]);
      InitLevelSetFromModel(fnl, center, extent, rotation, path, m);
    }
  } else if (name == "list") {
    if (sem()) {
      // TODO revise with bcast
      std::stringstream in;
      {
        std::stringstream path(var.String["eb_list_path"]);
        std::string filename;
        path >> filename;
        if (filename == "inline") {
          if (verb) {
            std::cerr << "InitLevelSet: Reading inline list of primitives from "
                         "list_path"
                      << std::endl;
          }
          in << path.rdbuf();
        } else {
          filename = path.str();
          std::ifstream fin(filename);
          if (verb) {
            std::cerr << "InitLevelSet: Open list of primitives '" << filename
                      << std::endl;
          }
          fassert(fin.good(), "Can't open list of primitives");
          in << fin.rdbuf();
        }
      }

      const size_t edim = var.Int["dim"];
      auto pp = UPrimList<Vect>::GetPrimitives(in, edim);

      if (verb) {
        std::cerr << "Read " << pp.size() << " primitives." << std::endl;
      }

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
    }
  } else {
    fassert(false, "unknown bc eb_init=" + name);
  }
  if (sem()) {
    if (var.Int["eb_init_inverse"]) {
      for (auto n : m.AllNodes()) {
        fnl[n] = -fnl[n];
      }
    }
    fnl.SetHalo(2);
  }
}

// Computes gradient at embeded face with Dirichlet conditions
// with second order using quadratic interpolation.
// rf: face center
// uf: given field value from boundary conditions
// nf: normal to face, towards excluded domain
// c: cell containign embedded face
// fcu: input field
template <class M, class T>
T GradDirichletQuadSecond(
    typename M::Vect rf, T uf, typename M::Vect nf, IdxCell c,
    const FieldCell<T>& fcu, const Embed<M>& eb) {
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
  auto& m = eb.GetMesh();
  auto h = m.GetCellSize()[0];
  auto dx = m.direction(nf.abs().argmax()).orient(-nf);
  auto dy = m.direction(1 - nf.abs().argmax());

  auto cw = m(c) + dx;
  auto dw = (rf[dx] - cw.center[dx]) / nf[dx];
  auto yw = rf[dy] - dw * nf[dy];

  auto cq = cw + dx;
  auto dq = (rf[dx] - cq.center[dx]) / nf[dx];
  auto yq = rf[dy] - dq * nf[dy];

  if (eb.IsExcluded(cw - dy) || eb.IsExcluded(cw) || eb.IsExcluded(cw + dy) ||
      eb.IsExcluded(cq - dy) || eb.IsExcluded(cq) || eb.IsExcluded(cq + dy)) {
    return GetNan<T>();
  }
  auto uw = interp::Quad(
      (yw - cw.center[dy]) / h, fcu[cw - dy], fcu[cw], fcu[cw + dy]);
  auto uq = interp::Quad(
      (yq - cq.center[dy]) / h, fcu[cq - dy], fcu[cq], fcu[cq + dy]);

  // XXX generated by src/gen/interp.py
  return (uf * sqr(dq) - uf * sqr(dw) + uq * sqr(dw) - uw * sqr(dq)) /
         (dq * dw * (dq - dw));
}

// Computes gradient at embeded face with Dirichlet conditions
// with first order using quadratic interpolation.
// rf: face center
// uf: given field value from boundary conditions
// nf: normal to face, towards excluded domain
// c: cell containign embedded face
// fcu: input field
template <class M, class T>
T GradDirichletQuad(
    typename M::Vect rf, T uf, typename M::Vect nf, IdxCell c,
    const FieldCell<T>& fcu, const Embed<M>& eb) {
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
  auto& m = eb.GetMesh();
  auto h = m.GetCellSize()[0];
  auto dx = m.direction(nf.abs().argmax()).orient(-nf);
  auto dy = m.direction(1 - nf.abs().argmax());
  auto cw = m(c) + dx;
  auto dist = (rf[dx] - cw.center[dx]) / nf[dx];
  auto yi = rf[dy] - dist * nf[dy];
  if (eb.IsExcluded(cw - dy) || eb.IsExcluded(cw) || eb.IsExcluded(cw + dy)) {
    return GetNan<T>();
  }
  auto ui = interp::Quad(
      (yi - cw.center[dy]) / h, fcu[cw - dy], fcu[cw], fcu[cw + dy]);
  return (uf - ui) / dist;
}

// Computes gradient at embeded face with Dirichlet conditions
// with first order using linear interpolation.
// rf: face center
// uf: given field value from boundary conditions
// nf: normal to face, towards excluded domain
// c: cell containign embedded face
// fcu: input field
template <class M, class T>
T GradDirichletLinear(
    typename M::Vect rf, T uf, typename M::Vect nf, IdxCell c,
    const FieldCell<T>& fcu, const Embed<M>& eb) {
  //            ---------------------            //
  //            |         |         |            //
  //            |         |  cw+dy  |            //
  //            |         |         |            //
  //            |\--------|---------|            //
  //            | \       |         |            //
  //            |  f  c   | cw=c+dx |            //
  //            |   \     |         |            //
  //            |----\----|---------|            //
  //                                             //
  auto& m = eb.GetMesh();
  auto h = m.GetCellSize()[0];
  auto dx = m.direction(nf.abs().argmax()).orient(-nf);
  auto dy = m.direction(1 - nf.abs().argmax()).orient(-nf);
  auto cw = m(c) + dx;
  auto dist = (rf[dx] - cw.center[dx]) / nf[dx];
  auto yi = rf[dy] - dist * nf[dy];
  if (eb.IsExcluded(cw) || eb.IsExcluded(cw + dy)) {
    return GetNan<T>();
  }
  auto ui = interp::Linear(
      dy.sign() * (yi - cw.center[dy]) / h, fcu[cw], fcu[cw + dy]);
  return (uf - ui) / dist;
}

// Computes gradient at embeded face with Dirichlet conditions
// with first order using linear fit.
// rf: face center
// uf: given field value from boundary conditions
// nf: normal to face, towards excluded domain
// c: cell containign embedded face
// fcu: input field
template <class M, class T>
T GradDirichletLinearFit(
    typename M::Vect rf, T uf, typename M::Vect nf, IdxCell c,
    const FieldCell<T>& fcu, const Embed<M>& eb) {
  auto h = eb.GetCellSize()[0];
  const T u = EvalLinearFit(rf - nf * h, c, fcu, eb);
  return (uf - u) / h;
}

template <class Vect_>
struct CombineMixed {
  using Vect = Vect_;
  using Scal = typename Vect::Scal;
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

  for (auto f : eb.SuFacesM()) {
    feu[f] = (fcu[f.cp] + fcu[f.cm]) * 0.5;
  }
  feu.GetFieldFace() = InterpolateBilinearFaces(feu.GetFieldFace(), eb);

  auto calc = [&](auto cf, IdxCell c, const BCond<T>& bc) {
    const Scal h = eb.GetCellSize()[0];
    switch (bc.type) {
      case BCondType::dirichlet: {
        return bc.val;
      }
      case BCondType::neumann: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const T u = EvalLinearFit(x, c, fcu, eb);
        return u + bc.val * (bc.nci == 0 ? 1 : -1) * h;
      }
      case BCondType::mixed:
      case BCondType::reflect: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const T u = EvalLinearFit(x, c, fcu, eb);
        return CombineMixed<Vect>()(bc.val, u + bc.val * h, eb.GetNormal(cf));
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

  for (auto f : eb.SuFacesM()) {
    feg[f] = (fcu[f.cp] - fcu[f.cm]) / h;
  }
  feg.GetFieldFace() = InterpolateBilinearFaces(feg.GetFieldFace(), eb);

  auto calc = [&](auto cf, IdxCell c, const BCond<T>& bc) {
    switch (bc.type) {
      case BCondType::dirichlet: {
        return (bc.nci == 0 ? 1 : -1) *
               GradDirichletLinearFit(
                   eb.GetFaceCenter(cf), bc.val, eb.GetNormal(cf), c, fcu, eb);
        //       GradDirichletQuadSecond(
        //       GradDirichletQuad(
        //       GradDirichletLinear(
      }
      case BCondType::neumann: {
        return bc.val * (bc.nci == 0 ? 1 : -1);
      }
      case BCondType::mixed:
      case BCondType::reflect: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const T u = EvalLinearFit(x, c, fcu, eb);
        return CombineMixed<Vect>()((bc.val - u) / h, bc.val, eb.GetNormal(cf));
      }
      case BCondType::extrap: {
        // TODO replace with dot-product of gradient and normal
        auto p = FitLinear(c, fcu, eb);
        const Vect x1 = eb.GetFaceCenter(cf);
        const Vect x0 = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const T u1 = ULinearFit<typename M::Vect>::EvalLinear(p, x1);
        const T u0 = ULinearFit<typename M::Vect>::EvalLinear(p, x0);
        return (u1 - u0) / h;
      }
    }
    return GetNan<T>();
  };
  mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto bc) { //
    feg[cf] = calc(cf, c, bc);
  });
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
        feu[f] = fcu[cm] + 0.5 * Superbee(du, -4 * fcg[cm].dot(rm) - du);
      } else if (fev[f] < 0) {
        const auto rp = m.GetVectToCell(f, 1);
        feu[f] = fcu[cp] - 0.5 * Superbee(du, 4 * fcg[cp].dot(rp) - du);
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
    switch (bc.type) {
      case BCondType::dirichlet: {
        return bc.val;
      }
      case BCondType::neumann: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const Scal u = EvalLinearFit(x, c, fcu, eb);
        return u + bc.val * (bc.nci == 0 ? 1 : -1) * h;
      }
      case BCondType::mixed:
      case BCondType::reflect: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const Scal u = EvalLinearFit(x, c, fcu, eb);
        return CombineMixed<Vect>()(bc.val, u + bc.val * h, eb.GetNormal(cf));
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
  for (auto f : eb.FacesM()) {
    const Scal sgn = (fev[f] > 0 ? 1 : -1);
    const auto c = (sgn > 0 ? f.cm() : f.cp()); // upwind cell
    auto d = f.direction();

    // temporal derivative, initial from acceleration
    Scal ut = (fc_src[f.cm] + fc_src[f.cp]) * 0.5;
    // normal derivative
    Scal ux = 0;
    // convective terms
    // normal direction
    {
      const IdxFace fm = c.face(-d);
      const IdxFace fp = c.face(d);
      // TODO: check why boundary faces treated differently in regular version
      ux += (is_boundary[fm] || is_boundary[fp]) ? feg[f]
                                                 : (feg[fm] + feg[fp]) * 0.5;
      const Scal w = fev[f] / eb.GetArea(f); // velocity
      ut -= ux * w;
    }
    // tangential directions
    for (size_t i = 0; i + 1 < M::dim; ++i) {
      d = d.next();
      const IdxFace fm = c.face(-d);
      const IdxFace fp = c.face(d);
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
    switch (bc.type) {
      case BCondType::dirichlet: {
        return bc.val;
      }
      case BCondType::neumann: {
        const Vect x = eb.GetFaceCenter(cf) - eb.GetNormal(cf) * h;
        const Scal u = EvalLinearFit(x, c, fcu, eb);
        return u + bc.val * (bc.nci == 0 ? 1 : -1) * h;
      }
      default:
        fassert(false, "unknown bc type");
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
  for (auto f : m.FacesM()) {
    const Scal sgn = (ffv[f] > 0 ? 1 : -1);
    const auto c = (sgn > 0 ? f.cm() : f.cp()); // upwind cell
    auto d = f.direction();

    // temporal derivative, initial from acceleration
    Scal ut = (fc_src[f.cm] + fc_src[f.cp]) * 0.5;
    // normal derivative
    Scal ux = 0;
    // convective terms
    // normal direction
    {
      const IdxFace fm = c.face(-d);
      const IdxFace fp = c.face(d);
      ux += (ffg[fm] + ffg[fp]) * 0.5;
      const Scal w = ffv[f] / m.GetArea(f); // velocity
      ut -= ux * w;
    }
    // tangential directions
    for (size_t i = 0; i + 1 < M::dim; ++i) {
      d = d.next();
      const IdxFace fm = c.face(-d);
      const IdxFace fp = c.face(d);
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
  auto calc = [&](IdxFace, IdxCell c, const BCond<Scal>& bc) {
    switch (bc.type) {
      case BCondType::dirichlet: {
        return bc.val;
      }
      case BCondType::neumann: {
        return fcu[c] + bc.val * (h * 0.5);
      }
      default:
        fassert(false, "unknown bc type");
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
  const Scal h = m.GetCellSize()[0];

  for (auto f : m.SuFacesM()) {
    ffu[f] = (fcu[f.cp] + fcu[f.cm]) * 0.5;
  }

  auto calc = [&](IdxFace f, IdxCell c, const BCond<T>& bc) {
    switch (bc.type) {
      case BCondType::dirichlet: {
        return bc.val;
      }
      case BCondType::neumann: {
        return fcu[c] + bc.val * (h * 0.5);
      }
      case BCondType::mixed:
      case BCondType::reflect: {
        const Scal q = (bc.nci == 0 ? 1. : -1.);
        const Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * q;
        return CombineMixed<Vect>()(
            bc.val, fcu[c] + bc.val * a, m.GetNormal(f));
      }
      case BCondType::extrap: {
        const IdxFace fo = m.GetFace(c, m.GetOpposite(m.GetNci(c, f)));
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
        fassert(false, "unknown bc type");
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

  for (auto f : m.SuFacesM()) {
    ffu[f] = (fcu[f.cp] - fcu[f.cm]) / h;
  }

  auto calc = [&](IdxFace f, IdxCell c, const BCond<T>& bc) {
    switch (bc.type) {
      case BCondType::dirichlet: {
        const Scal sgn = (bc.nci == 0 ? 1 : -1);
        const T uf = bc.val;
        const T u1 = fcu[c];
        const auto cnci = m.GetNci(c, f); // f=m.GetFace(c, cnci)
        const T u2 = fcu[m.GetCell(c, m.GetOpposite(cnci))];
        //    |      |      |
        //  uf|  u1  |  u2  |
        //    |      |      |
        return (uf * 8 - u1 * 9 + u2) / (3 * h) * sgn;
      }
      case BCondType::neumann: {
        return bc.val * (bc.nci == 0 ? 1 : -1);
      }
      case BCondType::extrap: {
        const Scal sgn = (bc.nci == 0 ? 1 : -1);
        const T u1 = fcu[c];
        const auto cnci = m.GetNci(c, f); // f=m.GetFace(c, cnci)
        const T u2 = fcu[m.GetCell(c, m.GetOpposite(cnci))];
        return (u1 - u2) / h * sgn;
      }
      default:
        fassert(false, "unknown bc type");
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
  const Scal h = m.GetCellSize()[0];

  if (scheme == ConvSc::superbee) {
    for (auto f : m.SuFaces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const Scal du = fcu[cp] - fcu[cm];
      if (ffv[f] > 0) {
        const auto rm = m.GetVectToCell(f, 0);
        ffu[f] = fcu[cm] + 0.5 * Superbee(du, -4 * fcg[cm].dot(rm) - du);
      } else if (ffv[f] < 0) {
        const auto rp = m.GetVectToCell(f, 1);
        ffu[f] = fcu[cp] - 0.5 * Superbee(du, 4 * fcg[cp].dot(rp) - du);
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

  auto calc = [&](IdxFace, IdxCell c, const BCond<Scal>& bc) -> Scal {
    switch (bc.type) {
      case BCondType::dirichlet: {
        return bc.val;
      }
      case BCondType::neumann: {
        return fcu[c] + bc.val * (h * 0.5);
      }
      default:
        fassert(false, "unknown bc type");
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
      const Scal wc = 1 / std::abs(eb.GetFaceOffset(c));
      T sum = feu[c] * wc;
      Scal sumw = wc;
      for (auto q : eb.Nci(c)) {
        IdxFace f = eb.GetFace(c, q);
        const Scal wf = 1 / std::abs(eb.GetFaceOffset(c, q));
        sum += feu[f] * wf;
        sumw += wf;
      }
      fcu[c] = sum / sumw;
    }
  }
  return fcu;
}

template <class M>
template <class T>
FieldCell<T> UEmbed<M>::Interpolate(const FieldFace<T>& feu, const M& m) {
  FieldCell<T> fcu(m, T(0));
  for (auto c : m.AllCells()) {
    T sum(0);
    Scal sumw = 0;
    for (auto q : m.Nci(c)) {
      sum += feu[m.GetFace(c, q)];
      sumw += 1;
    }
    fcu[c] = sum / sumw;
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
      auto p = ULinearFit<Vect>::FitLinear(xx, uu);
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
  auto excess = [&](IdxCell c) -> Scal {
    Scal dtmax = std::numeric_limits<Scal>::max();
    eb.LoopNci(c, [&](auto q) {
      const auto cf = eb.GetFace(c, q);
      const Scal flux = std::abs(ffv[cf]);
      if (flux != 0) {
        dtmax = std::min(dtmax, cfl * eb.GetVolume(c) / flux);
      }
    });
    return fcs[c] * std::max<Scal>(0, (dt - dtmax) / dt);
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
      auto dx = (dz.next(1)).orient(-n);
      auto dy = (dz.next(2)).orient(-n);
      auto r = (eb.GetFaceCenter(f) - f.center).abs();
      ffb[f] = interp::Bilinear(
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
  for (auto f : eb.FacesM()) {
    ExprFace e(0);
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
    e[2] -= fcu[f.cm] * e[0] + fcu[f.cp] * e[1];
    fee[f] = e;
  }

  auto calc = [&](auto cf, IdxCell c, const BCond<Scal>& bc) {
    ExprFace e(0);
    switch (bc.type) {
      case BCondType::dirichlet: {
        e[2] = feu[cf];
        break;
      }
      case BCondType::neumann: {
        e[bc.nci] = 1;
        e[2] = feu[cf];
        break;
      }
      default:
        fassert(false, "unknown bc type");
    }
    e[2] -= fcu[c] * e[bc.nci];
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
  const Scal h = m.GetCellSize()[0];

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

  auto calc = [&](IdxFace, IdxCell, const BCond<Scal>& bc) {
    ExprFace e(0);
    switch (bc.type) {
      case BCondType::dirichlet: {
        e[2] = bc.val;
        break;
      }
      case BCondType::neumann: {
        e[bc.nci] = 1;
        e[2] = bc.val * (h * 0.5);
        break;
      }
      default:
        fassert(false, "unknown bc type");
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
    const MapEmbed<BCond<Scal>>& mebc, const EB& eb) -> FieldEmbed<ExprFace> {
  FieldEmbed<ExprFace> feg(eb, ExprFace(0));
  const Scal h = eb.GetCellSize()[0];

  for (auto f : eb.Faces()) {
    const Scal a = 1 / h;
    feg[f] = ExprFace(-a, a, 0);
  }

  mebc.LoopBCond(eb, [&](auto cf, IdxCell, const auto& bc) {
    ExprFace e(0);
    switch (bc.type) {
      case BCondType::dirichlet: {
        const Scal a = 2 * (bc.nci == 0 ? 1 : -1) / h;
        e[bc.nci] = -a;
        e[2] = bc.val * a;
        break;
      }
      case BCondType::neumann: {
        e[2] = bc.val * (bc.nci == 0 ? 1 : -1);
        break;
      }
      default:
        fassert(false, "unknown bc type");
    }
    feg[cf] = e;
  });
  return feg;
}

template <class M>
auto UEmbed<M>::GradientImplicit(const MapEmbed<BCond<Scal>>& mebc, const M& m)
    -> FieldFace<ExprFace> {
  FieldFace<ExprFace> ffe(m, ExprFace(0));
  const Scal h = m.GetCellSize()[0];

  for (auto f : m.Faces()) {
    ExprFace e(0);
    const Scal a = 1 / h;
    ffe[f] = ExprFace(-a, a, 0);
  }

  // boundary conditions
  for (auto& p : mebc.GetMapFace()) {
    const auto f = p.first;
    const auto& bc = p.second;
    ExprFace e(0);
    switch (bc.type) {
      case BCondType::dirichlet: {
        const Scal a = 2 * (bc.nci == 0 ? 1 : -1) / h;
        e[bc.nci] = -a;
        e[2] = bc.val * a;
        break;
      }
      case BCondType::neumann: {
        e[2] = bc.val * (bc.nci == 0 ? 1 : -1);
        break;
      }
      default:
        fassert(false, "unknown bc type");
    }
    ffe[f] = e;
  }
  return ffe;
}

template <class M>
auto UEmbed<M>::GradientImplicit(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, const EB& eb)
    -> FieldEmbed<ExprFace> {
  const auto fee = Gradient(fcu, mebc, eb);
  auto feg = GradientImplicit(mebc, eb);
  eb.LoopFaces(
      [&](auto cf) { feg[cf].back() += fee[cf] - Eval(feg[cf], cf, fcu, eb); });
  return feg;
}

template <class M>
auto UEmbed<M>::GradientImplicit(
    const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, const M& m)
    -> FieldFace<ExprFace> {
  const auto ffe = Gradient(fcu, mebc, m);
  auto ffg = GradientImplicit(mebc, m);
  for (auto f : m.Faces()) {
    ffg[f].back() += ffe[f] - Eval(ffg[f], f, fcu, m);
  }
  return ffg;
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

template <class T, class MEB>
void Smoothen(
    FieldCell<T>& fc, const MapEmbed<BCond<T>>& mfc, MEB& eb, size_t iters) {
  auto& m = eb.GetMesh();
  auto sem = m.GetSem("smoothen");
  using UEB = UEmbed<typename MEB::M>;
  for (size_t i = 0; i < iters; ++i) {
    if (sem()) {
      fc = UEB::Interpolate(UEB::Interpolate(fc, mfc, eb), eb);
      m.Comm(&fc);
    }
    if (sem()) {
      // FIXME empty stage to finish communication in outer blocks
    }
  }
}

template <class M>
template <class MEB>
auto UEmbed<M>::GetVortScal(
    const FieldCell<Vect>& fcvel, const MapEmbed<BCond<Vect>>& me_vel,
    const MEB& eb) -> FieldCell<Scal> {
  auto& m = eb.GetMesh();

  std::array<FieldCell<Vect>, dim> grad;
  for (size_t d = 0; d < dim; ++d) {
    grad[d].Reinit(m, Vect(0));
    const auto mebc = GetScalarCond(me_vel, d, m);
    const FieldCell<Scal> fcu = GetComponent(fcvel, d);
    const FieldFace<Scal> ffg = Gradient(fcu, mebc, m);
    grad[d] = AverageGradient(ffg, m);
  }

  FieldCell<Scal> res(m, 0);
  for (auto c : m.Cells()) {
    res[c] = grad[1][c][0] - grad[0][c][1];
  }
  return res;
}

template <class M>
template <class MEB>
auto UEmbed<M>::GetVort(
    const FieldCell<Vect>& fcvel, const MapEmbed<BCond<Vect>>& me_vel,
    const MEB& eb) -> FieldCell<Vect> {
  auto& m = eb.GetMesh();

  FieldCell<Vect> r(m, Vect(0));
  if (M::dim < 3) {
    return r;
  }

  std::array<FieldCell<Vect>, M::dim> grad;
  for (size_t d = 0; d < M::dim; ++d) {
    grad[d].Reinit(m, Vect(0));
    const auto mebc = GetScalarCond(me_vel, d, m);
    const FieldCell<Scal> fcu = GetComponent(fcvel, d);
    const FieldFace<Scal> ffg = Gradient(fcu, mebc, m);
    grad[d] = AverageGradient(ffg, m);
  }

  for (auto c : eb.Cells()) {
    if (eb.IsExcluded(c)) continue;
    r[c][0] = grad[2][c][1] - grad[1][c][2];
    r[c][1] = grad[0][c][2] - grad[2][c][0];
    r[c][2] = grad[1][c][0] - grad[0][c][1];
  }
  return r;
}

// Converts vector conditions to scalar.
// mfv: vector velocity conditions
// d: direction, 0..2
template <class M>
MapEmbed<BCond<typename M::Scal>> GetScalarCond(
    const MapEmbed<BCond<typename M::Vect>>& mev, size_t d, const M& m) {
  using Scal = typename M::Scal;
  MapEmbed<BCond<Scal>> mes;

  for (auto& p : mev.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bcv = p.second;
    auto& bcs = mes[f];
    bcs.type = bcv.type;
    bcs.nci = bcv.nci;
    switch (bcv.type) {
      case BCondType::dirichlet:
      case BCondType::neumann:
        bcs.val = bcv.val[d];
        break;
      case BCondType::mixed:
      case BCondType::reflect:
        if (size_t(m.GetDir(f)) == d) {
          bcs.type = BCondType::dirichlet;
        } else {
          bcs.type = BCondType::neumann;
        }
        bcs.val = bcv.val[d];
        break;
      case BCondType::extrap:
        // nop
        break;
    }
  }
  for (auto& p : mev.GetMapCell()) {
    const IdxCell c = p.first;
    const auto& bcv = p.second;
    auto& bcs = mes[c];
    bcs.type = bcv.type;
    bcs.nci = bcv.nci;
    switch (bcv.type) {
      case BCondType::dirichlet:
      case BCondType::neumann:
        bcs.val = bcv.val[d];
        break;
      case BCondType::mixed:
      case BCondType::reflect:
        bcs.type = BCondType::neumann;
        // TODO revise to have zero normal component
      case BCondType::extrap:
        // nop
        break;
    }
  }
  return mes;
}
