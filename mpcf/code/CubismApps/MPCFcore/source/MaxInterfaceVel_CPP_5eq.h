/*
 *  MaxInterfaceVel_CPP_5eq.h
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer
 *  Date: October 2015
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#ifndef MAXINTERFACEVEL_CPP_5EQ_H
#define MAXINTERFACEVEL_CPP_5EQ_H

class MaxInterfaceVel_CPP_5eq
{
public:
    MaxInterfaceVel_CPP_5eq() {}

    virtual Real compute(const Real * const src, const int gptfloats) const;

    static void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NT, const int NBLOCKS, float MEASUREDTIME, const bool bAwk=false)
    {cout << "printflops() not yet implemented for MaxInterfaceVel_CPP_5eq" << endl; return;};
};

#endif
