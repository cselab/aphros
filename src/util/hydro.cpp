#include "hydro.ipp"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

template FieldCell<Scal> GetBcField(MapCondFaceFluid& mf, const M& m);

template FieldCell<Vect> GetVort(
    const FieldCell<Vect>& fcv, const MapCondFace& mf, M& m);

template void InitVel(FieldCell<Vect>& fcv, const Vars& var, const M& m);

template void GetFluidFaceCond(
    const Vars& var, const M& m, MapCondFaceFluid& mff,
    MapCondFaceAdvection<Scal>& mfa);

template void AppendBodyCond(
    const FieldCell<bool>& fc, std::string str, const M& m, Scal clear0,
    Scal clear1, Scal inletcl, Scal fill_vf,
    MapCell<std::shared_ptr<CondCellFluid>>* mcf, MapCondFaceFluid& mff,
    MapCondFaceAdvection<Scal>& mfa);

template void GetFluidCellCond(
    const Vars& var, M& m,
    MapCell<std::shared_ptr<CondCellFluid>>& mcvel,
    std::pair<Scal, int>& pdist);

template void DumpBcFaces(
    const MapCondFaceAdvection<Scal>& mfa, const MapCondFaceFluid& mff,
    std::string fn, M& m);
