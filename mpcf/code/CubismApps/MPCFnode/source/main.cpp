/*
 *  main.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <string>
#include <iostream>

#include "Test_SteadyState.h"
#include "FlowStep_LSRK3.h"
#include "Test_Cloud.h"
#include "Test_CloudAcoustic.h"

using namespace std;

Simulation * sim = NULL;

int main (int argc, const char ** argv)
{
	ArgumentParser parser(argc, argv);

  if (parser.exist("simConf"))    parser.readFile(parser("simConf").asString());

  sim = new Test_Cloud<Grid_t,FlowStep_LSRK3>(parser);

	sim->setup();

  return 0;
}


