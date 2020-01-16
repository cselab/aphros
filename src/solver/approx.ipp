// Created by Petr Karnakov on 29.12.2019
// Copyright 2019 ETH Zurich

#include "approx.h"

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
      throw std::runtime_error("GetCoeff: invalid ConvSc");
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

// Interpolates from cells to inner faces.
// T: value type (Scal or Vect)
// fc: field cell [s]
// fcgp: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M>
void Interpolate(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const MapCondFace& mfc, const FieldFace<T>& ffw, const M& m, ConvSc sc,
    typename M::Scal th, FieldFace<T>& ff) {
  InterpolateI(fc, fcgp, ffw, m, sc, th, ff);
  InterpolateB(fc, mfc, ff, m);
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
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    Scal a = m.GetArea(f) / m.GetVolume(cp);
    ff[f] = (fc[cp] - fc[cm]) * a;
  }
}

// Explicit gradient on boundary faces.
// fc: field [s]
// mfc: face conditions
// Output:
// ff: normal gradient [i]
template <class M, class T>
void GradientB(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m,
    FieldFace<T>& ff) {
  using Scal = typename M::Scal;

  for (const auto& it : mfc) {
    IdxFace f = it.GetIdx();
    const CondFace* cb = it.GetValue().Get(); // cond base
    if (auto cd = dynamic_cast<const CondFaceGrad<T>*>(cb)) {
      ff[f] = cd->GetGrad();
    } else if (auto cd = dynamic_cast<const CondFaceVal<T>*>(cb)) {
      size_t id = cd->GetNci();
      IdxCell c = m.GetCell(f, id);
      Scal g = (id == 0 ? 1. : -1.);
      Scal hr = m.GetArea(f) / m.GetVolume(c);
      Scal a = hr * 2 * g;
      ff[f] = (cd->GetValue() - fc[c]) * a;
    } else {
      throw std::runtime_error("GradientB: unknown cond");
    }
  }
}

// Explicit gradient on inner faces.
// fc: field [s]
// Output:
// ff: normal gradient [i]
template <class M, class T>
void Gradient(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m,
    FieldFace<T>& ff) {
  GradientI(fc, m, ff);
  GradientB(fc, mfc, m, ff);
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
  using Vect = GVect<Scal, 3>;
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
  using Vect = GVect<Scal, 3>;
  // v: value
  // n: normal to face
  static Scal Get(Scal v, const Vect& /*n*/) {
    return v;
  }
  static Vect Get(const Vect& v, const Vect& n) {
    return v - n * (2. * n.dot(v));
  }
};

// Interpolation to faces with defined conditions.
// fc: field cell [i]
// mfc: face cond
// Output:
// ff: values updated on faces defined in mfc
template <class T, class M>
void InterpolateB(
    const FieldCell<T>& fc, const MapCondFace& mfc, FieldFace<T>& ff,
    const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  for (const auto& it : mfc) {
    IdxFace f = it.GetIdx();
    const CondFace* cb = it.GetValue().Get(); // cond base
    size_t nci = cb->GetNci();
    if (auto cd = dynamic_cast<const CondFaceVal<T>*>(cb)) {
      ff[f] = cd->GetValue();
    } else if (auto cd = dynamic_cast<const CondFaceGrad<T>*>(cb)) {
      IdxCell c = m.GetCell(f, nci);
      Scal w = (nci == 0 ? 1. : -1.);
      Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * w;
      ff[f] = fc[c] + cd->GetGrad() * a;
    } else if (dynamic_cast<const CondFaceExtrap*>(cb)) {
      // TODO test
      IdxCell c = m.GetCell(f, nci);
      size_t q = m.GetNci(c, f);
      size_t qo = m.GetOpposite(q);
      IdxFace fo = m.GetFace(c, qo);
      Vect n = m.GetNormal(f);
      // cell
      const T& v0 = fc[c];
      Scal x0 = 0.;
      // opposite face
      const T& v1 = ff[fo];
      Scal x1 = n.dot(m.GetCenter(fo) - m.GetCenter(c));
      // target
      Scal xt = n.dot(m.GetCenter(f) - m.GetCenter(c));

      ff[f] = UExtrap(xt, x0, v0, x1, v1);
    } else if (dynamic_cast<const CondFaceReflect*>(cb)) {
      // TODO test
      IdxCell c = m.GetCell(f, nci);
      Vect n = m.GetNormal(f);
      auto v = fc[c];
      ff[f] = UReflectFace<Scal>::Get(v, n);
    } else {
      // TODO add name to CondFace etc
      throw std::runtime_error("InterpolateB: unknown cond");
    }
  }
}

// Interpolates from cells to support faces.
// T: value type (Scal or Vect)
// fc: field cell [a]
// mfc: face cond
// Output:
// field face [s]
template <class T, class M>
FieldFace<T> Interpolate(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m) {
  FieldFace<T> ff(m); // Valid 0 needed for CondFaceExtrap

  InterpolateS(fc, ff, m);
  InterpolateB(fc, mfc, ff, m);

  return ff;
}

template <class Scal>
Scal Superbee(Scal p, Scal q) {
  if (p > 0. && q > 0.) {
    return std::max(std::min(2 * p, q), std::min(p, 2 * q));
  } else if (p < 0. && q < 0.) {
    return -std::max(std::min(-2 * p, -q), std::min(-p, -2 * q));
  }
  return 0.;
}

// Second order upwind interpolation with TVD Superbee limiter
// fc: fieldcell [a]
// fcg: gradient of field [a]
// mfc: face cond
// ffw: flow direction [s]
// Output:
// fieldface [s]
template <class M>
FieldFace<typename M::Scal> InterpolateSuperbee(
    const FieldCell<typename M::Scal>& fc,
    const FieldCell<typename M::Vect>& fcg, const MapCondFace& mfc,
    const FieldFace<typename M::Scal>& ffw, const M& m, typename M::Scal th) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  FieldFace<Scal> ff(m);

  for (IdxFace f : m.SuFaces()) {
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    Vect rm = m.GetVectToCell(f, 0);
    Vect rp = m.GetVectToCell(f, 1);
    const auto& u = fc;
    const auto& g = fcg;
    Scal du = u[cp] - u[cm];
    if (ffw[f] > th) {
      ff[f] = u[cm] + 0.5 * Superbee(du, -4. * g[cm].dot(rm) - du);
    } else if (ffw[f] < -th) {
      ff[f] = u[cp] - 0.5 * Superbee(du, 4. * g[cp].dot(rp) - du);
    } else {
      ff[f] = 0.5 * (u[cm] + u[cp]);
    }
  }

  InterpolateB(fc, mfc, ff, m);
  return ff;
}

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

// Smoothens fieldcell.
// fc: fieldcell [s]
// mfc: condface
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void Smoothen(FieldCell<T>& fc, const MapCondFace& mfc, M& m, size_t rep) {
  auto sem = m.GetSem("smoothen");
  for (size_t i = 0; i < rep; ++i) {
    if (sem()) {
      fc = Average(Interpolate(fc, mfc, m), m);
      m.Comm(&fc);
    }
    // FIXME empty stage, without it cubismnc fails
    // on sim25 with m="128 16 16" np=2 OMP_NUM_THREADS=1 on two nodes
    // which is a minimal case with inner/halo blocks and MPI communication
    if (sem()) {
    }
  }
}

// Smoothens fieldcell with node-based averaging.
// fc: fieldcell [s]
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void SmoothenNode(FieldCell<T>& fc, M& m, size_t rep) {
  auto sem = m.GetSem("smoothen");
  for (size_t i = 0; i < rep; ++i) {
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
// computation: GetValue(field, idx) vs GetNorm(field)

template <class Field, class M, class Scal = typename M::Scal>
Scal CalcDiff(const Field& fa, const Field& fb, const M& m) {
  Scal r = 0.;
  using Idx = typename Field::Idx;
  for (Idx i : m.template GetIn<Idx>()) {
    r = std::max<Scal>(r, std::abs(fa[i] - fb[i]));
  }
  return r;
}

template <class Idx, class M, class Scal = typename M::Scal>
Scal CalcDiff(
    const GField<typename M::Vect, Idx>& fa,
    const GField<typename M::Vect, Idx>& fb, const M& m) {
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
std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z) {
  // TODO: test

  size_t s = z.size();
  std::vector<Scal> k(s);
  for (size_t i = 0; i < s; ++i) {
    Scal a = 0.;
    Scal b = 1.;
    for (size_t j = 0; j < s; ++j) {
      if (j != i) {
        b *= z[i] - z[j];
        Scal t = 1.;
        for (size_t k = 0; k < s; ++k) {
          if (k != i && k != j) {
            t *= x - z[k];
          }
        }
        a += t;
      }
    }
    k[i] = a / b;
  }
  return k;
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
void BcApply(FieldCell<T>& uc, const MapCondFace& mfc, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  for (const auto& it : mfc) {
    IdxFace f = it.GetIdx();
    auto& cb = it.GetValue();
    Vect n = m.GetNormal(f);
    IdxCell cmm, cm, cp, cpp;
    GetCellColumn(m, f, cb->GetNci(), cmm, cm, cp, cpp);
    if (cb.Get<CondFaceReflect>()) {
      uc[cm] = UReflectCell<Scal>::Get(uc[cp], n);
      uc[cmm] = UReflectCell<Scal>::Get(uc[cpp], n);
    } else if (auto cd = cb.Get<CondFaceVal<T>>()) {
      uc[cm] = cd->GetValue();
      uc[cmm] = cd->GetValue();
    }
  }
}

// Apply reflection on all boundaries
// fill: value for other types that CondFaceReflect
template <class T, class M>
void BcReflectAll(FieldCell<T>& uc, const MapCondFace& mfc, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  for (const auto& it : mfc) {
    IdxFace f = it.GetIdx();
    auto& cb = it.GetValue();
    Vect n = m.GetNormal(f);
    IdxCell cmm, cm, cp, cpp;
    GetCellColumn(m, f, cb->GetNci(), cmm, cm, cp, cpp);
    uc[cm] = UReflectCell<Scal>::Get(uc[cp], n);
    uc[cmm] = UReflectCell<Scal>::Get(uc[cpp], n);
  }
}
