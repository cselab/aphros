#include "hydro.ipp"


using M = MeshStructured<double, 3>;

template FieldCell<typename M::Scal> GetBcField(
    MapCondFaceFluid& mf, const M& m);

template FieldCell<typename M::Vect> GetVort(
    const FieldCell<typename M::Vect>& fcv, const MapCondFace& mf, M& m);

template void InitVel(FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m);

template void GetFluidFaceCond(
    const Vars& var, const M& m, MapCondFaceFluid& mfvel, MapCondFace& mfvf);

template void GetFluidCellCond(
    const Vars& var, M& m,
    MapCell<std::shared_ptr<solver::CondCellFluid>>& mcvel,
    std::pair<typename M::Scal, int>& pdist);

template void DumpBcFaces(
    const MapCondFace& mfc, const MapCondFaceFluid& mfcf, std::string fn, M& m);
