/*
 *  Types.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <fstream>
#include "math.h"
#include <iomanip>
#include <sstream>
#include <limits>
#include <vector>
#include <utility>
#include <iostream>
#include <omp.h>
#include <algorithm>

using namespace std;

#include <Grid.h>
#include <GridMorton.h>
#include <BlockLab.h>
#include <Profiler.h>
#include <ArgumentParser.h>
#include <SerializerIO.h>
#include <Timer.h>
#include <SerializerIO_VP.h>
#include <Timer.h>
#ifdef _USE_HDF_
#include <HDF5Dumper.h>
#endif
#include "GridTypes.h"

/* #define SETUP_MARKERS_IC \ */
/* const double mix_gamma = 1 + (G2*G1)/(G1*bubble+G2*(1-bubble)); \ */
/* const double mix_pinf  = (mix_gamma-1)/mix_gamma * (F1/G1*(1-bubble) + F2/G2*bubble); \ */
/* b(ix, iy, iz).G  = 1./(mix_gamma-1); \ */
/* b(ix, iy, iz).P = mix_gamma*mix_pinf/(mix_gamma-1); \ */
/* const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho; \ */
/* b(ix, iy, iz).energy   = pressure*b(ix, iy, iz).G + b(ix, iy, iz).P + ke; */

struct sort_pred {
    bool operator()(const std::pair<Real,Real> &left, const std::pair<Real,Real> &right) {
        return abs(left.first) < abs(right.first);
    }
};

class Simulation_Environment
{
public:
    static Real extent;
    static Real extents [3];
    static bool BC_PERIODIC [3];
};

class Simulation
{
public:
	virtual ~Simulation() {}
	virtual void run() = 0;
	virtual void setup() = 0;
};


struct StreamerGridPoint //dummy
{
    static const int channels = 1;

	void operate(const FluidElement& input, Real output[1]) const
	{
  		output[0] = input.alpha2;
	}
};

template <> inline void FluidBlock::Write<StreamerGridPoint>(
    ofstream& output, StreamerGridPoint streamer) const
{
    output.write((const char *)&data[0][0][0], 
        sizeof(FluidElement)*sizeX*sizeY*sizeZ);
}

template <> inline void FluidBlock::Read<StreamerGridPoint>(
    ifstream& input, StreamerGridPoint streamer)
{
    input.read((char *)&data[0][0][0], 
        sizeof(FluidElement)*sizeX*sizeY*sizeZ);
}
