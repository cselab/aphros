// Created by Petr Karnakov on 29.12.2019
// Copyright 2019 ETH Zurich

#include "approx.h"

// Retuns coefficients of interpolation scheme.
//   f = GetFace(cm, 1)
//   uf = u[cmm]*a[0] + u[cm]*a[1] + [cp]*a[2]
//    -------------------------------            //
//    |         |         |         |            //
//    |   cmm   |   cm    |f   cp   |            //
//    |         |         |         |            //
//    -------------------------------            //
// sc: scheme:
//   - fou: first order upwind
//   - cd: central differences (mean value)
//   - sou: second order upwind
//   - quick: QUICK
template <class Scal>
std::array<Scal, 3> GetCoeff(ConvSc sc) {
  std::array<Scal, 3> a;
  switch (sc) {
    case ConvSc::fou:
      a = {0., 1., 0.};
      break;
    case ConvSc::cd:
      a = {0., 0.5, 0.5};
      break;
    case ConvSc::sou:
      a = {-0.5, 1.5, 0.};
      break;
    case ConvSc::quick:
      a = {-1. / 8., 6. / 8., 3. / 8.};
      break;
    default:
      throw std::runtime_error(FILELINE + ": GetCoeff: invalid ConvSc");
  }
  return a;
}

// Explicit interpolation to faces near inner cells.
// fc: field cell [s]
// fc: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M>
void InterpolateI(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const FieldFace<T>& ffw, const M& m, ConvSc sc, typename M::Scal th,
    FieldFace<T>& ff) {
  using Scal = typename M::Scal;

  // f = fmm*a[0] + fm*a[1] + fp*a[2]
  std::array<Scal, 3> a = GetCoeff<Scal>(sc);

  ff.Reinit(m);
  for (auto f : m.Faces()) {
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    if (ffw[f] > th) {
      ff[f] = 4. * a[0] * fcgp[cm].dot(m.GetVectToCell(f, 0)) + a[1] * fc[cm] +
              (a[2] + a[0]) * fc[cp];
    } else if (ffw[f] < -th) {
      ff[f] = 4. * a[0] * fcgp[cp].dot(m.GetVectToCell(f, 1)) + a[1] * fc[cp] +
              (a[2] + a[0]) * fc[cm];
    } else {
      ff[f] = (fc[cm] + fc[cp]) * 0.5;
    }
  }
}

// Implicit interpolation to inner faces with deferred correction.
// fc: field cell [s]
// fc: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
//   - fou: first order upwind
//   - cd: central differences (mean value)
//   - sou: second order upwind
//   - quick: QUICK
// df: deferred correction factor (1. fully deferred)
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M, class Expr>
void InterpolateI(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const FieldFace<T>& ffw, FieldFace<Expr>& ff, const M& m, ConvSc sc,
    typename M::Scal df, typename M::Scal th) {
  using Scal = typename M::Scal;

  // f = fmm*a[0] + fm*a[1] + fp*a[2]
  std::array<Scal, 3> a = GetCoeff<Scal>(sc);

  ff.Reinit(m);
  auto dfm = 1. - df;
  for (auto f : m.Faces()) {
    Expr& e = ff[f];
    e.Clear();
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    if (ffw[f] > th) {
      e.InsertTerm(a[1] * dfm + df, cm);
      e.InsertTerm((a[2] + a[0]) * dfm, cp);
      e.SetConstant(
          4. * a[0] * fcgp[cm].dot(m.GetVectToCell(f, 0)) +
          (a[1] - 1.) * df * fc[cm] + (a[2] + a[0]) * df * fc[cp]);
    } else if (ffw[f] < -th) {
      e.InsertTerm((a[2] + a[0]) * dfm, cm);
      e.InsertTerm(a[1] * dfm + df, cp);
      e.SetConstant(
          4. * a[0] * fcgp[cp].dot(m.GetVectToCell(f, 1)) +
          (a[1] - 1.) * df * fc[cp] + (a[2] + a[0]) * df * fc[cm]);
    } else {
      e.InsertTerm(0.5, cm);
      e.InsertTerm(0.5, cp);
    }
  }
}

// Explicit gradient on faces near inner cells.
// fc: field [s]
// Output:
// ff: normal gradient [i]
template <class M, class T>
void GradientI(const FieldCell<T>& fc, const M& m, FieldFace<T>& ff) {
  using Scal = typename M::Scal;

  ff.Reinit(m);
  for (auto f : m.Faces()) {
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    const Scal a = m.GetArea(f) / m.GetVolume(cp);
    ff[f] = (fc[cp] - fc[cm]) * a;
  }
}

// Implicit gradient in inner faces.
// Output:
// ff: normal gradient [i]
template <class M, class Expr>
void GradientI(FieldFace<Expr>& ff, const M& m) {
  using Scal = typename M::Scal;

  ff.Reinit(m);
  for (auto f : m.Faces()) {
    Expr& e = ff[f];
    e.Clear();

    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    // XXX: adhoc uniform
    Scal a = m.GetArea(f) / m.GetVolume(cp);
    e.InsertTerm(-a, cm);
    e.InsertTerm(a, cp);
  }
}

// Interpolates from nodes to faces
template <class T, class M>
FieldFace<T> Interpolate(const FieldNode<T>& fn, const M& m) {
  FieldFace<T> ff(m);
  for (auto f : m.SuFaces()) {
    T s(0);
    for (size_t i = 0; i < m.GetNumNodes(f); ++i) {
      s += fn[m.GetNode(f, i)];
    }
    ff[f] = s / m.GetNumNodes(f);
  }
  return ff;
}

// Interpolation to inner faces.
// fc: field cell [s]
// Output:
// ff: face cell [i]
template <class T, class M>
void InterpolateI(const FieldCell<T>& fc, FieldFace<T>& ff, const M& m) {
  using Scal = typename M::Scal;

  for (auto f : m.Faces()) {
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    Scal a = 0.5;
    ff[f] = fc[cm] * (1. - a) + fc[cp] * a;
  }
}

// Interpolation to support faces.
// fc: field cell [a]
// Output:
// ff: face cell [s]
template <class T, class M>
void InterpolateS(const FieldCell<T>& fc, FieldFace<T>& ff, const M& m) {
  using Scal = typename M::Scal;

  for (auto f : m.SuFaces()) {
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    Scal a = 0.5;
    ff[f] = fc[cm] * (1. - a) + fc[cp] * a;
  }
}

template <class Scal>
class UReflectFace {
 public:
  using Vect = generic::Vect<Scal, 3>;
  // v: value
  // n: normal to face
  static Scal Get(Scal v, const Vect& /*n*/) {
    return v;
  }
  static Vect Get(const Vect& v, const Vect& n) {
    return v - n * n.dot(v);
  }
};

template <class Scal>
class UReflectCell {
 public:
  using Vect = generic::Vect<Scal, 3>;
  // v: value
  // n: normal to face
  static Scal Get(Scal v, const Vect& /*n*/) {
    return v;
  }
  static Vect Get(const Vect& v, const Vect& n) {
    return v - n * (2. * n.dot(v));
  }
};

// Returns average of fieldface.
// ff: fieldface [a]
// Output:
// fieldcell [a]
template <class T, class M>
FieldCell<T> Average(const FieldFace<T>& ff, const M& m) {
  using Scal = typename M::Scal;
  FieldCell<T> fc(m);
  for (IdxCell c : m.AllCells()) {
    T s(0);
    for (auto q : m.Nci(c)) {
      IdxFace f = m.GetFace(c, q);
      s += ff[f];
    }
    fc[c] = s / Scal(m.GetNumFaces(c));
  }
  return fc;
}

// Smoothens fieldcell with node-based averaging.
// fc: fieldcell [s]
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void SmoothenNode(FieldCell<T>& fc, M& m, size_t iters) {
  auto sem = m.GetSem("smoothen");
  for (size_t iter = 0; iter < iters; ++iter) {
    if (sem()) {
      using Scal = typename M::Scal;
      // generated by gen/smooth.py
      std::array<Scal, 27> a = {
          0.015625, 0.03125,  0.015625, 0.03125,  0.0625,   0.03125, 0.015625,
          0.03125,  0.015625, 0.03125,  0.0625,   0.03125,  0.0625,  0.125,
          0.0625,   0.03125,  0.0625,   0.03125,  0.015625, 0.03125, 0.015625,
          0.03125,  0.0625,   0.03125,  0.015625, 0.03125,  0.015625};

      using MIdx = typename M::MIdx;
      auto& bc = m.GetIndexCells();
      GBlock<IdxCell, M::dim> bo(MIdx(-1), MIdx(3));
      auto fcm = fc;
      for (auto c : m.Cells()) {
        MIdx w = bc.GetMIdx(c);
        T u = 0;
        size_t i = 0;
        for (MIdx wo : bo) {
          IdxCell cn = bc.GetIdx(w + wo);
          u += fcm[cn] * a[i++];
        }
        fc[c] = u;
      }
      m.Comm(&fc);
    }
  }
}

template <class T, class M>
void CommNodes(FieldNode<T>& fn, M& m) {
  auto sem = m.GetSem("commface");
  struct {
    std::array<FieldCell<T>, M::kCellNumNeighbourNodes> vfc;
  } * ctx(sem);
  auto& vfc = ctx->vfc;
  const auto range = GRange<size_t>(M::kCellNumNeighbourNodes);
  if (sem("comm")) {
    for (auto q : range) {
      vfc[q].Reinit(m);
    }
    for (auto c : m.Cells()) {
      for (auto q : range) {
        vfc[q][c] = fn[m.GetNode(c, q)];
      }
    }
    for (auto q : range) {
      m.Comm(&vfc[q]);
    }
  }
  if (sem("copy")) {
    for (auto c : m.AllCells()) {
      for (auto q : range) {
        fn[m.GetNode(c, q)] = vfc[q][c];
      }
    }
  }
  if (sem()) {
    // FIXME: empty stage required to prevent destruction of ctx
    // until communication is finished in outer blocks
  }
}

// Smoothens fieldcell with node-based averaging.
// fc: fieldcell [s]
// iters: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void SmoothenNode(FieldNode<T>& fn, M& m, size_t iters) {
  auto sem = m.GetSem("smoothen-node");
  for (size_t iter = 0; iter < iters; ++iter) {
    if (sem()) {
      using Scal = typename M::Scal;
      // generated by gen/smooth.py
      std::array<Scal, 27> weight = {
          0.015625, 0.03125,  0.015625, 0.03125,  0.0625,   0.03125, 0.015625,
          0.03125,  0.015625, 0.03125,  0.0625,   0.03125,  0.0625,  0.125,
          0.0625,   0.03125,  0.0625,   0.03125,  0.015625, 0.03125, 0.015625,
          0.03125,  0.0625,   0.03125,  0.015625, 0.03125,  0.015625};

      using MIdx = typename M::MIdx;
      auto& index = m.GetIndexNodes();
      GBlock<IdxCell, M::dim> stencil(MIdx(-1), MIdx(3));
      auto fnm = fn;
      for (auto n : m.Nodes()) {
        const MIdx w = index.GetMIdx(n);
        T u = 0;
        size_t i = 0;
        for (MIdx wo : stencil) {
          const IdxNode nn = index.GetIdx(w + wo);
          u += fnm[nn] * weight[i++];
        }
        fn[n] = u;
      }
    }
    if (sem.Nested()) {
      CommNodes(fn, m);
    }
  }
}

// Returns gradient.
// ff: scalar fieldface [s]
// Output:
// gradient vector [s]
template <class M>
FieldCell<typename M::Vect> Gradient(
    const FieldFace<typename M::Scal>& ff, const M& m) {
  using Vect = typename M::Vect;
  FieldCell<Vect> fc(m);
  for (auto c : m.SuCells()) {
    Vect s(0);
    for (auto q : m.Nci(c)) {
      IdxFace f = m.GetFace(c, q);
      s += m.GetOutwardSurface(c, q) * ff[f];
    }
    fc[c] = s / m.GetVolume(c);
  }
  return fc;
}

// Convention: Use Get/Set for fast procedures and Calc for those requiring
// computation: second(field, idx) vs GetNorm(field)

template <class Field, class M>
typename M::Scal CalcDiff(const Field& fa, const Field& fb, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0.;
  using Idx = typename Field::Idx;
  for (Idx i : m.template GetIn<Idx>()) {
    r = std::max<Scal>(r, std::abs(fa[i] - fb[i]));
  }
  return r;
}

template <class Idx, class M>
typename M::Scal CalcDiff(
    const GField<typename M::Vect, Idx>& fa,
    const GField<typename M::Vect, Idx>& fb, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0.;
  for (Idx i : m.template GetIn<Idx>()) {
    r = std::max<Scal>(r, (fa[i] - fb[i]).norminf());
  }
  return r;
}

// Coefficients for approximation of gradient with polynomial.
// x: target point
// z: stencil points
// Output:
// k: such that grad(x) = sum_i (ki * f(zi))
template <class Scal>
std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& stencil) {
  // TODO: test

  size_t size = stencil.size();
  std::vector<Scal> res(size);
  for (size_t i = 0; i < size; ++i) {
    Scal a = 0;
    Scal b = 1;
    for (size_t j = 0; j < size; ++j) {
      if (j != i) {
        b *= stencil[i] - stencil[j];
        Scal t = 1;
        for (size_t k = 0; k < size; ++k) {
          if (k != i && k != j) {
            t *= x - stencil[k];
          }
        }
        a += t;
      }
    }
    res[i] = a / b;
  }
  return res;
}

// Returns GetGradCoeffs(x,z[b:]) preceeded by b zeros.
template <class Scal>
std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z, size_t b) {
  size_t s = z.size();
  size_t ss = s - b;
  std::vector<Scal> zz(ss);
  for (size_t i = 0; i < ss; ++i) {
    zz[i] = z[b + i];
  }
  std::vector<Scal> kk = GetGradCoeffs(x, zz);
  std::vector<Scal> k(s);
  for (size_t i = 0; i < b; ++i) {
    k[i] = 0.;
  }
  for (size_t i = 0; i < ss; ++i) {
    k[b + i] = kk[i];
  }
  return k;
}

// Apply boudnary conditions to halo cells
template <class T, class M>
void BcApply(FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  for (const auto& p : me.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    const Vect n = m.GetNormal(f);
    IdxCell cmm, cm, cp, cpp;
    GetCellColumn(m, f, bc.nci, cmm, cm, cp, cpp);
    if (bc.type == BCondType::reflect) {
      uc[cm] = UReflectCell<Scal>::Get(uc[cp], n);
      uc[cmm] = UReflectCell<Scal>::Get(uc[cpp], n);
    } else if (bc.type == BCondType::dirichlet) {
      uc[cm] = bc.val;
      uc[cmm] = bc.val;
    } else if (bc.type == BCondType::extrap) {
    }
  }
}

// Apply reflection on all boundaries
// fill: value for other types that CondFaceReflect
template <class T, class M>
void BcReflectAll(
    FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  for (const auto& p : me.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    const Vect n = m.GetNormal(f);
    IdxCell cmm, cm, cp, cpp;
    GetCellColumn(m, f, bc.nci, cmm, cm, cp, cpp);
    uc[cm] = UReflectCell<Scal>::Get(uc[cp], n);
    uc[cmm] = UReflectCell<Scal>::Get(uc[cpp], n);
  }
}
