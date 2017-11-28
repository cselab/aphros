/*
 *  Profiler.cpp
 *  Cubism
 *
 *  Created by Diego Rossinelli on 9/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include <sys/time.h>

#include "Profiler.h"

//#include <tbb/tick_count.h>
//using namespace tbb;

void ProfileAgent::_getTime(ClockTime& time)
{
	//time = tick_count::now();
	gettimeofday(&time, NULL);
}

double ProfileAgent::_getElapsedTime(const ClockTime& tS, const ClockTime& tE)
{
	return (tE.tv_sec - tS.tv_sec) + 1e-6 * (tE.tv_usec - tS.tv_usec);
	//return (tE - tS).seconds();
}



