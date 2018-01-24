/*
 *  FlowStep_LSRK3.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include "FlowStep_LSRK3.h"

#include <cstdlib>
#include <omp.h>

using namespace std;

namespace LSRK3data {
    int verbosity;
    int step_id = 0;
}


