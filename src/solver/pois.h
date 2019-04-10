#pragma once

#include <exception>
#include <memory>
#include <cmath>
#include <sstream>
#include <string>

#include "geom/mesh.h"
#include "linear/linear.h"
#include "cond.h"
#include "solver/approx.h"


namespace solver {

// Solves Poisson equation: \nabla \nabla u = r.
// fcr: rhs [i]
// fce: buffer for system
// mf: boundary conditions, replaced with zero-gradient
// Output:
// fcu: solution [a]
template <class M>
class PoisSolver {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;


  PoisSolver(const MapFace<std::shared_ptr<solver::CondFace>>& mf, M& m) 
      : m(m), mf_(mf) {}
  // Solve linear system fce = 0
  // fce: expressions [i]
  // Output:
  // fc: solution [a]
  // m.GetSolveTmp(): modified temporary fields
  void LinSolve(const FieldCell<Expr>& fce, FieldCell<Scal>& fc, M& m) {
    auto sem = m.GetSem("solve");
    if (sem("solve")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);
      auto l = ConvertLs(fce, *lsa, *lsb, *lsx, m);
      using T = typename M::LS::T; 
      l.t = T::symm; // solver type
      l.prefix = "vort"; // XXX: adhoc for vorticity
      m.Solve(l);
    }
    if (sem("copy")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);

      fc.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc[c] = (*lsx)[i++];
      }
      CHECKNAN(fc, m.CN());
      m.Comm(&fc);
    }
  }
  void Solve(const FieldCell<Scal>& fcr) {
    auto sem = m.GetSem("pois");
    auto& fcrt = fcu_;      // temporary rhs
    if (sem("reduce")) {
      fcrt = fcr;             
      sumr_ = 0.;
      sumv_ = 0.;
      for (auto c : m.Cells()) {
        auto v = m.GetVolume(c);
        sumr_ += fcrt[c] * v;
        sumv_ += v;
      }
      m.Reduce(&sumr_, "sum");
      m.Reduce(&sumv_, "sum");
    }
    if (sem("assemble")) {
      // compute average rhs
      sumr_ /= sumv_;

      FieldFace<Expr> ffe(m); // normal derivative
      GradientI(ffe, m);
      // overwrite boundaries
      FaceGradB<M, Expr> gb(m, mf_);
      for (auto it : mf_) {
        IdxFace f = it.GetIdx();
        Expr e = gb.GetExpr(f);
        e.SortTerms();
      }

      fce_.Reinit(m);    // equations in cells
      for (auto c : m.Cells()) {
        auto& e = fce_[c];
        e.Clear();
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          e += ffe[f] * (m.GetOutwardFactor(c, q) * m.GetArea(f));
        }
        e += Expr(-(fcr[c] - sumr_) * m.GetVolume(c));
      }
    }
    if (sem.Nested("solve")) {
      LinSolve(fce_, fcu_, m);
    }
    if (sem("comm")) {
      m.Comm(&fcu_);
      fce_.Free();
    }
  }
  const FieldCell<Scal>& GetField() const {
    return fcu_;
  }

 private:
  M& m;
  FieldCell<Expr> fce_;
  FieldCell<Scal> fcu_;
  MapFace<std::shared_ptr<solver::CondFace>> mf_;
  Scal sumr_; // sum of rhs * volume
  Scal sumv_; // sum of volume
};


} // namespace solver

