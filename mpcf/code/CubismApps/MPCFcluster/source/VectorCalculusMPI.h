/*
 *  VectorCalculusMPI.h
 *  MPCFcluster
 *
 *  Created by Fabian Wermelinger 06/24/2015
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */
#ifndef VECTORCALCULUSMPI_H_LXKGZA6P
#define VECTORCALCULUSMPI_H_LXKGZA6P

#include "FlowStep_LSRK3MPI.h"

// generic operator evaluator
template <typename TLab, typename TKernel, typename TGrid>
inline void evaluateOperator_MPI(TGrid& grid)
{
    TKernel kappa;
    process<TLab>(kappa, grid, 0, 0);
};

// miscellaneous operators
// ----------------------------------------------------------------------------
#include "VectorOperator.h"

#endif /* VECTORCALCULUSMPI_H_LXKGZA6P */
