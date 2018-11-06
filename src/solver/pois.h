#pragma once

#include <exception>
#include <memory>
#include <cmath>
#include <sstream>
#include <string>

#include "geom/mesh.h"
#include "linear/linear.h"
#include "cond.h"


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


  PoisSolver(const MapFace<std::shared_ptr<solver::CondFace>>& mf, M& m) {
    // zero-derivative bc for Scal
    for (auto it : mf_velcond_) {
      IdxFace i = it.GetIdx();
      mfz_[i] = std::make_shared<
          CondFaceGradFixed<Scal>>(0, it.GetValue()->GetNci());
    }
  }
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
  void Solve(FieldCell<T>& fcr) {
    auto sem = m.GetSem("pois");
    if (sem("assemble")) {
      const FieldFace<Expr> ffe(m); // normal derivative
      // set all faces
      for (auto f : m.Faces()) {
        auto& e = ffe[f];
        e.Clear();
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        Vect dm = m.GetVectToCell(f, 0);
        Vect dp = m.GetVectToCell(f, 1);
        Scal a = m.GetArea / (dp - dm).norm();
        e.InsertTerm(-a, cm);
        e.InsertTerm(a, cp);
      }
      // overwrite boundaries
      for (auto it : mf) {
        IdxFace f = it.GetIdx();
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        auto& e = ffe[f];
        e.Clear();
        e.InsertTerm(0, cm);
        e.InsertTerm(0, cp);
      }

      fce_.Reinit(m);
      for (auto c : m.Cells()) {
        auto& e = fce_[c];
        e.Clear();
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          e += ffe[f] * m.GetOutwardFactor(c, q);
        }
        e += Expr(fcr[c] * m.GetVolume(c));
      }
    }
    if (sem.Nested("solve")) {
      LinSolve(fce_, fcu_, m)
    }
  }

 private:
  FieldCell<Expr>& fce_;
  FieldCell<Scal>& fcu_;
  MapFace<std::shared_ptr<solver::CondFace>> mfz_;
};


} // namespace solver

