#include "fluid.h"
#include "fluid.ipp"


using M = MeshStructured<double, 3>;

template FieldCell<typename M::Scal> GetBcField(
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mf, const M& m);

template FieldCell<typename M::Vect> GetVort(
    const FieldCell<typename M::Vect>& fcv,
    const MapFace<std::shared_ptr<solver::CondFace>>& mf, M& m);

template void InitVel(FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m);

template void GetFluidFaceCond(
    const Vars& var, const M& m,
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mfvel,
    MapFace<std::shared_ptr<solver::CondFace>>& mfvf);

template void GetFluidCellCond(
    const Vars& var, M& m,
    MapCell<std::shared_ptr<solver::CondCellFluid>>& mcvel,
    std::pair<typename M::Scal, int>& pdist);
