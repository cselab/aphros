/*
 *  Tests.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/20/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "TestLabs.h"

//#elif defined(_BCLABCLOUD1DCHARNREF_ACOUSTIC_)
// 1d characteristic-based non-refelcting boundary condition with acoustic
// forcing pressure
#define _CHARACTERISTIC_1D_BOUNDARY_
typedef BlockLabCloud1DCharNonReflectAcousticForcing_5eq<Block_t, std::allocator> Lab;
