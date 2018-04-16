#pragma once

#include <exception>
#include <fstream>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"

namespace solver {

using namespace geom;

template <class M>
FieldCell<typename M::Vect> GetDeformingVelocity(const M& m) {
  using Vect = typename M::Vect;
  FieldCell<Vect> r(m, 0);
  for (auto c : m.Cells()) {
    auto x = m.GetCenter(c);
    r[c][0] = -std::cos(x[1]) * std::sin(x[0]);
    r[c][1] = std::cos(x[0]) * std::sin(x[1]);
  }
  return r;
}

template <class M>
class AdvectionSolver : public UnsteadyIterativeSolver {
  using Mesh = M;
  using Scal = typename Mesh::Scal;

 protected:
  Mesh& m;
  const FieldFace<Scal>* ffv_; // volume flux [velocity*area]
  const FieldCell<Scal>* fcs_; // source [value/time]

 public:
  AdvectionSolver(double t, double dt,
                  Mesh& m,
                  const FieldFace<Scal>* ffv /*volume flux*/,
                  const FieldCell<Scal>* fcs /*source*/)
      : UnsteadyIterativeSolver(t, dt)
      , m(m)
      , ffv_(ffv)
      , fcs_(fcs) {}
  virtual const FieldCell<Scal>& GetField(Layers) = 0;
  virtual const FieldCell<Scal>& GetField() {
    return GetField(Layers::time_curr);
  }
};


} // namespace solver
