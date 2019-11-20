#include "convdiffvg.ipp"
#include "convdiffi.h"
#include "convdiffe.h"

namespace solver {

template class ConvDiffVectGeneric<MeshStructured<double, 3>,
                               ConvDiffScalImp<MeshStructured<double, 3>>>;

template class ConvDiffVectGeneric<MeshStructured<double, 3>,
                               ConvDiffScalExp<MeshStructured<double, 3>>>;

} // namespace solver

