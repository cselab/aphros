#include <exception>
#include <memory>
#include <array>
#include <stdexcept>
#include <string>

namespace solver {

template <class Idx, class Expr>
class Approx {
 public:
  virtual ~Approx() {}
  virtual Expr GetExpr(Idx) const = 0;
};

// Convection scheme
enum class ConvSc { fou, cd, sou, quick };

std::string GetName(ConvSc sc) {
  switch (sc) {
    case ConvSc::fou: return "fou";
    case ConvSc::cd: return "cd";
    case ConvSc::sou: return "sou";
    case ConvSc::quick: return "quick";
    default: throw std::runtime_error("InterpolateI: invalid ConvSc");
  }
  return "";
}

ConvSc GetConvSc(std::string s) {
  auto l = {ConvSc::fou, ConvSc::cd, ConvSc::sou, ConvSc::quick};
  for (ConvSc sc : l) {
    if (GetName(sc) == s) {
      return sc;
    }
  }
  throw std::runtime_error("ConvSc: invalid name=" + s);
}

// Interpolation to inner faces with deferred correction.
// in terms of exprssions.
// fc: field cell [s]
// fc: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
//   - fou: first order upwind
//   - cd: central differences (mean value)
//   - sou: second order upwind
//   - quick: QUICK
// Output:
// ff: face cell [i]
template <class T, class M, class Expr>
void InterpolateI(const FieldCell<T>& fc,
                  const FieldCell<typename M::Vect>& fcgp,
                  const FieldFace<T>& ffw, FieldFace<Expr>& ff, 
                  const M& m, ConvSc sc, 
                  typename M::Scal th=1e-10) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // f = fmm*a[0] + fm*a[1] + fp*a[2]
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
      a = {-1./8., 6./8., 3./8.}; 
      break;
    default:
      throw std::runtime_error("InterpolateI: invalid ConvSc");
  }

  ff.Reinit(m);
  for (auto f : m.Faces()) {
    Expr& e = ff[f];
    e.Clear();
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    if (ffw[f] > th) {
      e.InsertTerm(1., cm);
      e.InsertTerm(0., cp);
      e.SetConstant(4. * a[0] * fcgp[cm].dot(m.GetVectToCell(f, 0)) +
          (a[1] - 1.) * fc[cm] + (a[2] + a[0]) * fc[cp]);
    } else if (ffw[f] < -th) {
      e.InsertTerm(0., cm);
      e.InsertTerm(1., cp);
      e.SetConstant(4. * a[0] * fcgp[cp].dot(m.GetVectToCell(f, 1)) +
          (a[1] - 1.) * fc[cp] + (a[2] + a[0]) * fc[cm]);
    } else {
      e.InsertTerm(0.5, cm);
      e.InsertTerm(0.5, cp);
    }
  }
}

template <class M, class Expr>
class FaceVal : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  FaceVal(const M& m, const FieldFace<Scal>& w, 
          const FieldCell<Scal>& fcp, const FieldCell<Vect>& fcgp, 
          Scal th = 1e-10)
      : m(m), w_(w), fcp_(fcp), fcgp_(fcgp), th_(th) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    // f = fmm + fm + fp
    //const std::array<Scal, 3> a = {0., 1., 0.}; // FOU
    //const std::array<Scal, 3> a = {0., 0.5, 0.5}; // CD
    //const std::array<Scal, 3> a = {-0.5, 1.5, 0.}; // SOU
    const std::array<Scal, 3> a = {-1./8., 6./8., 3./8.}; // QUICK
    if (w_[f] > th_) {
      e.InsertTerm(a[1], cm);
      e.InsertTerm(a[2] + a[0], cp);
      e.SetConstant(4. * a[0] * fcgp_[cm].dot(m.GetVectToCell(f, 0)));
    } else if (w_[f] < -th_) {
      e.InsertTerm(a[2] + a[0], cm);
      e.InsertTerm(a[1], cp);
      e.SetConstant(4. * a[0] * fcgp_[cp].dot(m.GetVectToCell(f, 1)));
    } else {
      e.InsertTerm(0.5, cm);
      e.InsertTerm(0.5, cp);
    }
    return e;
  }

 private:
  const M& m;
  const FieldFace<Scal>& w_; // flow direction (wind)
  const FieldCell<Scal>& fcp_; // previous
  const FieldCell<Vect>& fcgp_; // gradient of previous
  const Scal th_;
};

template <class M, class Expr>
class FaceValB : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  FaceValB(const M& m, const MapFace<std::shared_ptr<CondFace>>& mfc)
      : m(m), mfc_(mfc) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    if (auto cb = mfc_.find(f)) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      e.InsertTerm(0, cm);
      e.InsertTerm(0, cp);
      if (auto cd = dynamic_cast<CondFaceVal<Scal>*>(cb->get())) {
        e.SetConstant(cd->GetValue());
      } else if (auto cd = dynamic_cast<CondFaceGrad<Scal>*>(cb->get())) {
        size_t id = cd->GetNci();
        IdxCell cc = m.GetNeighbourCell(f, id);
        Scal g = (id == 0 ? 1. : -1.);
        Scal a = m.GetVectToCell(f, id).norm() * g;
        e.SetConstant(a * cd->GetGrad());
        e.InsertTerm(1., cc);
      } else {
        throw std::runtime_error("Unknown boundary condition type");
      }
    } else {
      throw std::runtime_error("Boundary condition not set");
    }
    return e;
  }

 private:
  const M& m;
  const MapFace<std::shared_ptr<CondFace>>& mfc_;
};

template <class M, class Expr>
class FaceGrad : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  explicit FaceGrad(const M& m) : m(m) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    auto dm = m.GetVectToCell(f, 0);
    auto dp = m.GetVectToCell(f, 1);
    Scal a = Scal(1) / (dp - dm).norm();
    e.InsertTerm(-a, cm);
    e.InsertTerm(a, cp);
    return e;
  }

 private:
  const M& m;
};

// Gradient to inner faces.
// fc: field cell [s]
// fc: gradient [s]
// ffw: flow direction [i]
// Output:
// ff: face cell [i]
template <class M, class Expr>
void GradientI(FieldFace<Expr>& ff, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  ff.Reinit(m);
  for (auto f : m.Faces()) {
    Expr& e = ff[f];
    e.Clear();

    IdxCell cm = m.GetNeighbourCell(f, 0);
    IdxCell cp = m.GetNeighbourCell(f, 1);
    auto dm = m.GetVectToCell(f, 0);
    auto dp = m.GetVectToCell(f, 1);
    Scal a = Scal(1) / (dp - dm).norm();
    e.InsertTerm(-a, cm);
    e.InsertTerm(a, cp);
  }
}


template <class M, class Expr>
class FaceGradB : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  explicit FaceGradB(
      const M& m, const MapFace<std::shared_ptr<CondFace>>& mfc)
      : m(m), mfc_(mfc) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    if (auto cb = mfc_.find(f)) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      e.InsertTerm(0, cm);
      e.InsertTerm(0, cp);
      if (auto cd = dynamic_cast<CondFaceGrad<Scal>*>(cb->get())) {
        e.SetConstant(cd->GetGrad());
      } else if (auto cd = dynamic_cast<CondFaceVal<Scal>*>(cb->get())) {
        size_t id = cd->GetNci();
        IdxCell c = m.GetNeighbourCell(f, id);
        Scal g = (id == 0 ? 1. : -1.);
        Scal a = 1. / m.GetVectToCell(f, id).norm() * g;
        e.SetConstant(a * cd->GetValue());
        e.InsertTerm(-a, c);
      } else {
        throw std::runtime_error("Unknown boundary condition type");
      }
    } else {
      throw std::runtime_error("Boundary condition not set");
    }
    return e;
  }

 private:
  const M& m;
  const MapFace<std::shared_ptr<CondFace>>& mfc_;
};

} // namespace solver

