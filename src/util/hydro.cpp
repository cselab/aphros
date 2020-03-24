// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "hydro.ipp"

using M = MeshStructured<double, 3>;
using EB = Embed<M>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

template FieldCell<Scal> GetBcField(MapCondFaceFluid& mf, const M& m);

template FieldCell<Vect> GetVort(
    const FieldCell<Vect>& fcv, const MapCondFace& mf, M& m);

template void InitVel(FieldCell<Vect>& fcv, const Vars& var, const M& m);

template void GetFluidFaceCond(
    const Vars& var, const M& m, MapCondFaceFluid& mff,
    MapCondFaceAdvection<Scal>& mfa);

template std::tuple<
    MapEmbed<BCondFluid<typename M::Vect>>,
    MapEmbed<CondFaceAdvection<typename M::Scal>>, MapEmbed<size_t>,
    std::vector<std::string>>
InitBc(const Vars& var, const M& eb);

template std::tuple<
    MapEmbed<BCondFluid<typename M::Vect>>,
    MapEmbed<CondFaceAdvection<typename M::Scal>>, MapEmbed<size_t>,
    std::vector<std::string>>
InitBc(const Vars& var, const EB& eb);

template void AppendBodyCond(
    const FieldCell<bool>& fc, std::string str, const M& m, Scal clear0,
    Scal clear1, Scal inletcl, Scal fill_vf,
    MapCell<std::shared_ptr<CondCellFluid>>* mcf, MapCondFaceFluid& mff,
    MapCondFaceAdvection<Scal>& mfa);

template void GetFluidCellCond(
    const Vars& var, M& m, MapCell<std::shared_ptr<CondCellFluid>>& mcvel);

template void DumpBcFaces(
    const MapCondFaceAdvection<Scal>& mfa, const MapCondFaceFluid& mff,
    std::string fn, M& m);

template void InitVort(
    const FieldCell<Vect>& fcvort, FieldCell<Vect>& fcvel,
    const MapCondFaceFluid& mf_fluid, bool verb, M& m);

template void DumpTraj(
    M& m, bool dm, const Vars& var, size_t frame, typename M::Scal t,
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
    const FieldCell<typename M::Scal>& fc_sig, const MapCondFace& mf_sig,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldCell<typename M::Scal>& fcvf,
    const FieldFace<typename M::Scal>& ffvfsm, const AdvectionSolver<M>* asb);

template void ProjectVolumeFlux(
    FieldFace<typename M::Scal>& ffv, const MapCondFaceFluid& mfc, M& m);

template std::map<Scal, Scal> CalcArea(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Vect>*> fcn,
    const Multi<const FieldCell<typename M::Scal>*> fca,
    const Multi<const FieldCell<typename M::Scal>*> fccl, M& m);

template std::map<Scal, Scal> CalcVolume(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*> fcu,
    const Multi<const FieldCell<typename M::Scal>*> fccl, M& m);
