// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "hydro.ipp"

using M = MeshStructured<double, 3>;
using EB = Embed<M>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

template FieldCell<Vect> GetVort(
    const FieldCell<Vect>& fcv, const MapEmbed<BCond<Vect>>& mebc, M& m);

template void InitVel(FieldCell<Vect>& fcv, const Vars& var, const M& m);

template std::tuple<
    MapEmbed<BCondFluid<typename M::Vect>>,
    MapEmbed<BCondAdvection<typename M::Scal>>, MapEmbed<size_t>,
    std::vector<std::string>,
    std::vector<std::map<std::string, typename M::Scal>>>
InitBc(const Vars& var, const M& eb, std::set<std::string> known_keys);

template std::tuple<
    MapEmbed<BCondFluid<typename M::Vect>>,
    MapEmbed<BCondAdvection<typename M::Scal>>, MapEmbed<size_t>,
    std::vector<std::string>,
    std::vector<std::map<std::string, typename M::Scal>>>
InitBc(const Vars& var, const EB& eb, std::set<std::string> known_keys);

template void GetFluidCellCond(
    const Vars& var, M& m, MapCell<std::shared_ptr<CondCellFluid>>& mcvel);

template void InitVort(
    const FieldCell<Vect>& fcvort, FieldCell<Vect>& fcvel,
    const MapEmbed<BCondFluid<Vect>>& mebc_fluid, M& m);

template void DumpTraj(
    M& m, bool dm, const Vars& var, size_t frame, typename M::Scal t,
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*>& fcvf,
    const Multi<const FieldCell<typename M::Scal>*>& fccl,
    const Multi<const FieldCell<typename M::MIdx>*>& fcim,
    const FieldCell<typename M::Scal>& fcp,
    const FieldCell<typename M::Vect>& fcvel,
    const FieldCell<typename M::Vect>& fcvelm, typename M::Scal dt);
template void DumpTraj(
    EB& m, bool dm, const Vars& var, size_t frame, typename M::Scal t,
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*>& fcvf,
    const Multi<const FieldCell<typename M::Scal>*>& fccl,
    const Multi<const FieldCell<typename M::MIdx>*>& fcim,
    const FieldCell<typename M::Scal>& fcp,
    const FieldCell<typename M::Vect>& fcvel,
    const FieldCell<typename M::Vect>& fcvelm, typename M::Scal dt);

template void CalcSurfaceTension(
    const M& m, const GRange<size_t>& layers, const Vars& var,
    FieldCell<typename M::Vect>& fc_force,
    FieldFace<typename M::Scal>& ff_force,
    const FieldCell<typename M::Scal>& fc_sig,
    const MapEmbed<BCond<typename M::Scal>>& mf_sig,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldCell<typename M::Scal>& fcvf,
    const FieldFace<typename M::Scal>& ffvfsm, const AdvectionSolver<M>* asb);

template void ProjectVolumeFlux(
    FieldFace<typename M::Scal>& ffv,
    const MapEmbed<BCondFluid<typename M::Vect>>& mfc, M& m);

template std::map<Scal, Scal> CalcArea(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Vect>*> fcn,
    const Multi<const FieldCell<typename M::Scal>*> fca,
    const Multi<const FieldCell<typename M::Scal>*> fccl, M& m);

template std::map<Scal, Scal> CalcVolume(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*> fcu,
    const Multi<const FieldCell<typename M::Scal>*> fccl, M& m);
