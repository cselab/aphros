#pragma once

#include "solver/fluid.h"
#include "parse/vars.h"

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
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mf, const M& m);

// Computes vorticity of vector field.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
FieldCell<typename M::Vect> GetVort(const FieldCell<typename M::Vect>& fcv,
                       const MapFace<std::shared_ptr<solver::CondFace>>& mf,
                       M& m);

// Initializes velocity from parameters.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
void InitVel(FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m);

// Returns fluid conditions on domain boundaries.
// Output:
// mfvel: conditions for velocity and pressure
// mfvf: conditions for volume fraction
template <class M>
void GetFluidFaceCond(
    const Vars& var, const M& m,
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mfvel,
    MapFace<std::shared_ptr<solver::CondFace>>& mfvf);


// Returns fluid cell conditions.
// Output:
// mcvel: output
// pdist, pdistmin: temporary buffer for reduction, 
// TODO: revise, allow temporary buffers in functions (attached to m)
template <class M>
void GetFluidCellCond(
    const Vars& var, M& m,
    MapCell<std::shared_ptr<solver::CondCellFluid>>& mcvel,
    std::pair<typename M::Scal, int>& pdist);
