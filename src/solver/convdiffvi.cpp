#include "convdiffvi.h"
#include "convdiffvi.ipp"
#include "convdiffi.h"
#include "convdiffe.h"

namespace solver {

template class ConvDiffVectImp<MeshStructured<double, 3>,
                               ConvDiffScalImp<MeshStructured<double, 3>>>;

template class ConvDiffVectImp<MeshStructured<double, 3>,
                               ConvDiffScalExp<MeshStructured<double, 3>>>;

} // namespace solver

