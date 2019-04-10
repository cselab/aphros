#pragma once

#include <memory>

#include "solver/cond.h"
#include "geom/map.h"


// Converts vector conditions to scalar.
// mfv: vector velocity conditions
// d: direction, 0..2
template <class M>
MapFace<std::shared_ptr<solver::CondFace>> 
GetScalarCond(const MapFace<std::shared_ptr<solver::CondFace>>& mfv, 
              size_t d, const M& m) {
  using namespace solver;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  MapFace<std::shared_ptr<CondFace>> mfs;
  // Face conditions for each velocity component
  for (auto it : mfv) {
    IdxFace f = it.GetIdx();
    CondFace* cb = it.GetValue().get();
    if (auto p = dynamic_cast<CondFaceVal<Vect>*>(cb)) {
      mfs[f] = std::make_shared<CondFaceValComp<Vect>>(p, d);
    } else if (auto p = dynamic_cast<CondFaceGrad<Vect>*>(cb)) {
      mfs[f] = std::make_shared<CondFaceGradComp<Vect>>(p, d);
    } else if (auto p = dynamic_cast<CondFaceReflect*>(cb)) {
      auto nci = cb->GetNci();
      // XXX: adhoc for cartesian grid
      if (d == m.GetNormal(f).abs().argmax()) { 
        // normal, zero value
        mfs[f] = std::make_shared<CondFaceValFixed<Scal>>(0., nci);
      } else { 
        // tangential, zero gradient
        mfs[f] = std::make_shared<CondFaceGradFixed<Scal>>(0., nci);
      }
    } else {
      throw std::runtime_error("GetScalarCond: unknown face condition");
    }
  }
  return mfs;
}
