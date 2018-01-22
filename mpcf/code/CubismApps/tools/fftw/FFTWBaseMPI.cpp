/*
 *  FFTWBaseMPI.cpp
 *
 *  Created by Ursula Rasthofer on 6/29/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#include "FFTWBaseMPI.h"

int FFTWBaseMPI::registered_objects_ = 0;
bool FFTWBaseMPI::initialized_ = false;
