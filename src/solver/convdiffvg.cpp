#include "convdiffvg.ipp"
#include "convdiffe.h"
#include "convdiffi.h"

template class ConvDiffVectGeneric<
    MeshStructured<double, 3>, ConvDiffScalImp<MeshStructured<double, 3>>>;

template class ConvDiffVectGeneric<
    MeshStructured<double, 3>, ConvDiffScalExp<MeshStructured<double, 3>>>;
