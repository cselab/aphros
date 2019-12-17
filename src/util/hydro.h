#pragma once

#include "parse/vars.h"
#include "solver/advection.h"
#include "solver/convdiff.h"
#include "solver/convdiffv.h"
#include "solver/fluid.h"

using namespace solver;

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
FieldCell<typename M::Scal> GetBcField(MapCondFaceFluid& mf, const M& m);

// Computes vorticity of vector field.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
FieldCell<typename M::Vect> GetVort(
    const FieldCell<typename M::Vect>& fcv, const MapCondFace& mf, M& m);

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
    const Vars& var, const M& m, MapCondFaceFluid& mff,
    MapCondFaceAdvection<typename M::Scal>& mfa);

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

// Appends step-wise approximation of body to cell and face conditions.
// Shape is defined as fc=1.
// Boundary conditions added on faces separating fc=0 and fc=1.
// Neighbor cell index (Nci) is chosen from cell fc=0.
// fc: boolean mask, fc=1 is inside the body
// Output:
// mcf: fluid cell conditions
// bc: boundary condition
// mff,mfa: fluid and advection face conditions
// pdist, pdistmin: temporary buffer for reduction,
// TODO: revise, allow temporary buffers in functions (attached to m)
template <class M>
void AppendBodyCond(
    const FieldCell<bool>& fc, std::string, M& m,
    MapCell<std::shared_ptr<solver::CondCellFluid>>& mcf,
    MapCondFaceFluid& mff, MapCondFaceAdvection<typename M::Scal>& mfa);

// Dumps faces with boundary conditions.
// mfc: boundary conditions
// fn: filename
template <class M>
void DumpBcFaces(
    const MapCondFaceAdvection<typename M::Scal>& mfa,
    const MapCondFaceFluid& mff, std::string fn, M& m);
