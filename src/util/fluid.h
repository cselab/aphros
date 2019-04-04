#pragma once

#include "solver/fluid.h"

// Returns field with the type (index)
// of boundary conditions in an adjacent face:
//   0: empty
//   1: no-slip wall
//   2: free-slip wall
//   3: inlet
//   4: outlet
//   -1: unknown
// mf: boundary conditions
template <class M>
FieldCell<typename M::Scal> GetBcField(
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mf, const M& m) {
  FieldCell<typename M::Scal> fc(m, 0);
  for (auto it : mf) {
    IdxFace f = it.GetIdx();
    auto* b = it.GetValue().get();
    size_t nci = b->GetNci();
    IdxCell c = m.GetNeighbourCell(f, nci);
    if (dynamic_cast<solver::fluid_condition::NoSlipWall<M>*>(b)) {
      fc[c] = 1.;
    } else if (dynamic_cast<solver::fluid_condition::SlipWall<M>*>(b)) {
      fc[c] = 2.;
    } else if (dynamic_cast<solver::fluid_condition::Inlet<M>*>(b)) {
      fc[c] = 3.;
    } else if (dynamic_cast<solver::fluid_condition::Outlet<M>*>(b)) {
      fc[c] = 4.;
    } else {
      fc[c] = -1.;
    }
  }
  return fc;
}

// Computes vorticity of vector field.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
FieldCell<typename M::Vect> GetVort(const FieldCell<typename M::Vect>& fcv,
                       const MapFace<std::shared_ptr<solver::CondFace>>& mf,
                       M& m) {
  auto ffv = solver::Interpolate(fcv, mf, m);

  auto d0 = solver::Gradient(GetComponent(ffv, 0), m);
  auto d1 = solver::Gradient(GetComponent(ffv, 1), m);
  auto d2 = solver::Gradient(GetComponent(ffv, 2), m);

  FieldCell<typename M::Vect> r(m);
  for (auto c : m.Cells()) {
    r[c][0] = d2[c][1] - d1[c][2];
    r[c][1] = d0[c][2] - d2[c][0];
    r[c][2] = d1[c][0] - d0[c][1];
  }

  return r;
}

