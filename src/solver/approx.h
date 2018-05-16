#include <exception>
#include <memory>
#include <array>

namespace solver {

template <class Idx, class Expr>
class Approx {
 public:
  virtual ~Approx() {}
  virtual Expr GetExpr(Idx) const = 0;
};

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
      } else if (auto cd =
          dynamic_cast<CondFaceGrad<Scal>*>(cb->get())) {
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
      if (auto cd = dynamic_cast<
          CondFaceGrad<Scal>*>(cb->get())) {
        e.SetConstant(cd->GetGrad());
      } else if (auto cd = dynamic_cast<
          CondFaceVal<Scal>*>(cb->get())) {
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

