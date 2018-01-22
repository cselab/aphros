/*
 *  VectorCalculus.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 10/13/2016
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */
#ifndef VECTORCALCULUSMPI_H_LXKGZA6P
#define VECTORCALCULUSMPI_H_LXKGZA6P

#include "BlockProcessor_OMP.h"

// generic operator evaluator
template <typename TLab, typename TKernel, typename TGrid>
inline void evaluateOperator_OMP(TGrid& grid)
{
    TKernel kappa;
    process_OMP<TLab>(kappa, grid, 0, 0);
};

// miscellaneous operators
// ----------------------------------------------------------------------------
#include "VectorOperator.h"

#endif /* VECTORCALCULUSMPI_H_LXKGZA6P */
