// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <map>
#include <set>

#include "parse/vars.h"
#include "solver/advection.h"
#include "solver/convdiff.h"
#include "solver/convdiffv.h"
#include "solver/fluid.h"
#include "solver/multi.h"

// Computes vorticity of vector field.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
FieldCell<typename M::Vect> GetVort(
    const FieldCell<typename M::Vect>& fcv,
    const MapEmbed<BCond<typename M::Vect>>& mebc, M& m);

// Initializes velocity from parameters.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
void InitVel(FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m);

// Reads boundary conditions. Error is generated for unknown keys.
// known_keys: other keys to be parsed and returned in map<string, Scal>
// Output:
// - fluid conditions BCondFluid on faces,
// - advection conditions BCondAdvection on faces,
// - group index on faces
// - string descriptors of boundary conditions for each group
// - custom keys and values allowed by known_keys for each group
template <class MEB>
std::tuple<
    MapEmbed<BCondFluid<typename MEB::Vect>>,
    MapEmbed<BCondAdvection<typename MEB::Scal>>, MapEmbed<size_t>,
    std::vector<std::string>,
    std::vector<std::map<std::string, typename MEB::Scal>>>
InitBc(const Vars& var, const MEB& eb, std::set<std::string> known_keys);

template <class MEB>
void DumpBcPoly(
    const std::string filename, const MapEmbed<size_t>& me_group,
    const MapEmbed<typename MEB::Scal>& me_contang, const MEB& meb,
    typename MEB::M& m);

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
    MapCell<std::shared_ptr<CondCellFluid>>* mcf,
    MapEmbed<BCondFluid<typename M::Vect>>& mff,
    MapEmbed<BCondAdvection<Scal>>& mfa);

// Computes velocity fcvel from vorticity fcvort
template <class M>
void InitVort(
    const FieldCell<typename M::Vect>& fcvort,
    FieldCell<typename M::Vect>& fcvel,
    const MapEmbed<BCondFluid<typename M::Vect>>& mebc_fluid,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m);

template <
    class EB, class Scal = typename EB::Scal, class Vect = typename EB::Vect>
void DumpTraj(
    EB& m, bool dm, const Vars& var, size_t frame, Scal t,
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcvf,
    const Multi<const FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<typename EB::MIdx>*>& fcim,
    const FieldCell<Scal>& fcp, const FieldCell<Vect>& fcvel,
    const FieldCell<Vect>& fcvelm, Scal dt);

template <class EB>
void CalcTraj(
    EB& eb, const GRange<size_t>& layers,
    const Multi<const FieldCell<typename EB::Scal>*>& fcvf,
    const Multi<const FieldCell<typename EB::Scal>*>& fccl,
    const Multi<const FieldCell<typename EB::MIdx>*>& fcim,
    const FieldCell<typename EB::Scal>& fcp,
    const FieldCell<typename EB::Vect>& fcvel,
    /*out*/
    std::vector<std::string>& column_names,
    std::vector<typename EB::Scal>& row_colors,
    std::vector<std::vector<typename EB::Scal>>& table);

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
    const FieldCell<typename M::Scal>& fc_sig,
    const MapEmbed<BCond<typename M::Scal>>& mf_sig,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldCell<typename M::Scal>& fcvf,
    const FieldFace<typename M::Scal>& ffvfsm, const AdvectionSolver<M>* as);

// ffv: volume flux
// mfc: face fluid conditions
// Output:
// ffv: divergence-free volume flux
template <class M>
void ProjectVolumeFlux(
    FieldFace<typename M::Scal>& ffv,
    const MapEmbed<BCondFluid<typename M::Vect>>& mfc,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m);

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
