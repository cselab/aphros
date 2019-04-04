#include "fluid.h"
#include "fluid.ipp"


using M = MeshStructured<double, 3>;

template FieldCell<typename M::Scal> GetBcField(
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mf, const M& m);

template FieldCell<typename M::Vect> GetVort(
    const FieldCell<typename M::Vect>& fcv,
    const MapFace<std::shared_ptr<solver::CondFace>>& mf, M& m);

template void InitVel(FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m);
