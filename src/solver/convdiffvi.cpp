#include "convdiffvi.h"
#include "convdiffvi.ipp"

namespace solver {

template class ConvDiffVectImp<MeshStructured<double, 3>>;

} // namespace solver

