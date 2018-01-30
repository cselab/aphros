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
//#include <GridMorton.h>
#include <BlockLab.h>
#include <Profiler.h>
#include <ArgumentParser.h>
#include <SerializerIO.h>
#include <Timer.h>
#include <HDF5Dumper.h>
#include "GridTypes.h"


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
