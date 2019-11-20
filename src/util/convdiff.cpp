#include "convdiff.h"

#include "solver/convdiffvg.h"
#include "solver/convdiffi.h"
#include "solver/convdiffe.h"

template <class M>
std::unique_ptr<ConvDiffVect<M>> GetConvDiff<M>::operator()(
    Conv conv, M& m, const FieldCell<Vect>& fcw,
    const MapFace<std::shared_ptr<CondFace>>& mfc,
    const MapCell<std::shared_ptr<CondCell>>& mcc,
    const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
    const FieldCell<Vect>* fcs, const FieldFace<Scal>* ffv,
    double t, double dt, const Par& par) {
  using CDI = ConvDiffVectGeneric<M, ConvDiffScalImp<M>>; // implicit
  using CDE = ConvDiffVectGeneric<M, ConvDiffScalExp<M>>; // explicit

  switch (conv) {
    case Conv::imp:
        return std::unique_ptr<CDI>(
            new CDI(m, fcw, mfc, mcc, fcr, ffd, fcs, ffv, t, dt, par));
    case Conv::exp:
        return std::unique_ptr<CDE>(
            new CDE(m, fcw, mfc, mcc, fcr, ffd, fcs, ffv, t, dt, par));
  }
  return nullptr;
}

using M = MeshStructured<double, 3>;

template std::unique_ptr<ConvDiffVect<M>> GetConvDiff<M>::operator()(
    Conv conv, M& m, const FieldCell<Vect>& fcw,
    const MapFace<std::shared_ptr<CondFace>>& mfc,
    const MapCell<std::shared_ptr<CondCell>>& mcc,
    const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
    const FieldCell<Vect>* fcs, const FieldFace<Scal>* ffv,
    double t, double dt, const Par& par);

