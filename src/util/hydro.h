// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include "parse/vars.h"
#include "solver/advection.h"
#include "solver/convdiff.h"
#include "solver/convdiffv.h"
#include "solver/fluid.h"
#include "solver/multi.h"

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

template <class MEB>
std::tuple<
    MapEmbed<BCondFluid<typename MEB::Vect>>,
    MapEmbed<CondFaceAdvection<typename MEB::Scal>>, MapEmbed<size_t>,
    std::vector<std::string>>
InitBc(const Vars& var, const MEB& eb);

// Returns fluid cell conditions.
// Output:
// mcvel: output
// pdist, pdistmin: temporary buffer for reduction,
// TODO: revise, allow temporary buffers in functions (attached to m)
template <class M>
void GetFluidCellCond(
    const Vars& var, M& m, MapCell<std::shared_ptr<CondCellFluid>>& mcvel);

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
template <class M, class Scal = typename M::Scal>
void AppendBodyCond(
    const FieldCell<bool>& fc, std::string str, const M& m, Scal clear0,
    Scal clear1, Scal inletcl, Scal fill_vf,
    MapCell<std::shared_ptr<CondCellFluid>>* mcf, MapCondFaceFluid& mff,
    MapCondFaceAdvection<Scal>& mfa);

// Dumps faces with boundary conditions.
// mfc: boundary conditions
// fn: filename
template <class M>
void DumpBcFaces(
    const MapCondFaceAdvection<typename M::Scal>& mfa,
    const MapCondFaceFluid& mff, std::string fn, M& m);

// Computes velocity fcvel from vorticity fcvort
template <class M, class Vect = typename M::Vect>
void InitVort(
    const FieldCell<Vect>& fcvort, FieldCell<Vect>& fcvel,
    const MapCondFaceFluid& mf_fluid, bool verb, M& m);

template <class M>
void DumpTraj(
    M& m, bool dm, const Vars& var, size_t frame, typename M::Scal t,
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*>& fcvf,
    const Multi<const FieldCell<typename M::Scal>*>& fccl,
    const Multi<const FieldCell<typename M::MIdx>*>& fcim,
    const FieldCell<typename M::Scal>& fcp,
    const FieldCell<typename M::Vect>& fcvel,
    const FieldCell<typename M::Vect>& fcvelm, typename M::Scal dt);

template <class M>
void AppendSurfaceTension(
    const M& m, FieldFace<typename M::Scal>& ffst,
    const FieldCell<typename M::Scal>& fcu,
    const FieldCell<typename M::Scal>& fck,
    const FieldFace<typename M::Scal>& ffsig);

template <class M>
void AppendSurfaceTension(
    const M& m, FieldFace<typename M::Scal>& ffst, const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*> fcu,
    const Multi<const FieldCell<typename M::Scal>*> fccl,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldFace<typename M::Scal>& ffsig);

// fc_force: force field to append
// ff_force: face force field to append
// fc_sig: surface tension coefficient
// mf_sig: boundary conditions for fc_sig
// fck: curvature
// fcvf: volume fraction
// ffvfsm: smoothed volume fraction
template <class M>
void CalcSurfaceTension(
    const M& m, const GRange<size_t>& layers, const Vars& var,
    FieldCell<typename M::Vect>& fc_force,
    FieldFace<typename M::Scal>& ff_force,
    const FieldCell<typename M::Scal>& fc_sig, const MapCondFace& mf_sig,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldCell<typename M::Scal>& fcvf,
    const FieldFace<typename M::Scal>& ffvfsm, const AdvectionSolver<M>* as);

// ffv: volume flux
// mfc: face fluid conditions
// Output:
// ffv: divergence-free volume flux
template <class M>
void ProjectVolumeFlux(
    FieldFace<typename M::Scal>& ffv, const MapCondFaceFluid& mfc, M& m);

// Returns total volume of each color, map[cl]=area
template <class M>
std::map<typename M::Scal, typename M::Scal> CalcArea(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Vect>*> fcn,
    const Multi<const FieldCell<typename M::Scal>*> fca,
    const Multi<const FieldCell<typename M::Scal>*> fccl, M& m);

// Returns total area of each color, map[cl]=volume
template <class M>
std::map<typename M::Scal, typename M::Scal> CalcVolume(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*> fcu,
    const Multi<const FieldCell<typename M::Scal>*> fccl, M& m);
