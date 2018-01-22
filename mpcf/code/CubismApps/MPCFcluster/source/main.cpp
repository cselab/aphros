/*
 *  main.cpp
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <mpi.h>
#ifdef _QPXEMU_
#include <xmmintrin.h>
#endif

#include <ArgumentParser.h>
#include <Timer.h>
//#include <BlockProcessing.h>

#include "Tests.h"
#include "Test_SteadyStateMPI.h"
#include "FlowStep_LSRK3MPI.h"

typedef FlowStep_LSRK3MPI<GridMPI_t> TFlowStep;

//#include "Test_ShockBubbleMPI.h"
//#include "Test_SICMPI.h"
#include "Test_CloudMPI.h"
#include "Test_CloudAcousticMPI.h"
#include "Test_SICCloudMPI.h"
// #include "Test_2DSBIMPI.h"
#include "Test_ChannelFlowMPI.h"
// #include "Test_HITMPI.h"
// //#include "Test_HITCloudMPI.h"
// #include "Test_TGVMPI.h"
// //#include "Test_KHinstabilityMPI.h"
// #include "Test_OscillatingShapeMPI.h"
#ifdef _USE_FFTW_
// #include "Test_CVTMPI.h"
// #include "Test_RingCloudMPI.h"
// #include "Test_HITMPI.h"
#endif

#ifdef _QPX_
#include <spi/include/l1p/sprefetch.h>
#endif

using namespace std;

Simulation * sim = NULL;

int main (int argc, const char ** argv)
{
//	{
//	char *ptr;
//	ptr = getenv("MPI_INIT_THREAD");
//	if (ptr != NULL) {
		int provided;
		MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &provided);	// peh: dropped C++ bindings
//	} else {
//		MPI_Init(&argc, (char ***)&argv);
//	}
//	}

	//L1P_SetStreamPolicy(L1P_stream_optimistic);

#ifdef _QPX_
	//peh: pick prefetching policy
        #pragma omp parallel
        {
	  int hwthread = Kernel_PhysicalHWThreadID(); // returns 0-3 on each core
	  char *ptr;
	  int stream_depth;

	  if (hwthread == 0) {
	    ptr = getenv("L1P_POLICY");
	    if (ptr != NULL) {
	      if (strncasecmp(ptr,"opt", 3) == 0) L1P_SetStreamPolicy(L1P_stream_optimistic);
	      else if (strncasecmp(ptr,"con", 3) == 0) L1P_SetStreamPolicy(L1P_stream_confirmed);
	      else if (strncasecmp(ptr,"dcbt",4) == 0) L1P_SetStreamPolicy(L1P_confirmed_or_dcbt);
	      else if (strncasecmp(ptr,"dis", 3) == 0) L1P_SetStreamPolicy(L1P_stream_disable);
	      else L1P_SetStreamPolicy(L1P_confirmed_or_dcbt);
	    }
	    else {
	      L1P_SetStreamPolicy(L1P_confirmed_or_dcbt);
	    }

	    ptr = getenv("L1P_DEPTH");
	    if (ptr != NULL) {
	      stream_depth = atoi(ptr);
	      L1P_SetStreamDepth(stream_depth);
	    }
	  }
        }
#endif

        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        const bool isroot = (0 == myrank);

#ifdef _GIT_HASH_
    const string git_hash(_GIT_HASH_);
#else
    const string git_hash("UNKNOWN");
#endif

	if (isroot)
    {
        cout << "GIT COMMIT: " << git_hash << endl;
		cout << "=================  MPCF cluster =================" << endl;
    }

	ArgumentParser parser(argc, argv);

    // read configuration files, if available
    if (parser.exist("cubismConf"))
    {
        parser.readFile(parser("cubismConf").asString());
        const size_t NX = parser("NX").asInt();
        const size_t NY = parser("NY").asInt();
        const size_t NZ = parser("NZ").asInt();
        const size_t PX = parser("xpesize").asInt();
        const size_t PY = parser("ypesize").asInt();
        const size_t PZ = parser("zpesize").asInt();
        const size_t BX = parser("bpdx").asInt();
        const size_t BY = parser("bpdy").asInt();
        const size_t BZ = parser("bpdz").asInt();
        const bool Xgood = (NX == _BLOCKSIZE_*BX*PX);
        const bool Ygood = (NY == _BLOCKSIZE_*BY*PY);
        const bool Zgood = (NZ == _BLOCKSIZE_*BZ*PZ);
        const bool notGood = !(Xgood && Ygood && Zgood);
        if (notGood)
        {
            if (isroot) cout << "Sanity check failed... Check file: " << parser("cubismConf").asString() << endl;
            abort();
        }
    }
    if (parser.exist("simConf"))    parser.readFile(parser("simConf").asString());

	const bool bFlush2Zero = parser("-f2z").asBool(true);

#ifdef _QPXEMU_
	if (bFlush2Zero)
#pragma omp parallel
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	}
#endif

	parser.unset_strict_mode();

        if (!isroot)
		parser.mute();

        if (isroot)
	  cout << "Dispatcher: " << parser("-dispatcher").asString("omp") << endl;

	parser.set_strict_mode();

	//Environment::setup(max(1, parser("-nthreads").asInt()));

    if( parser("-sim").asString() == "steady" )
        sim = new Test_SteadyStateMPI<GridMPI_t,TFlowStep>(MPI_COMM_WORLD, parser);
//    else if( parser("-sim").asString() == "sb" )
//        sim = new Test_ShockBubbleMPI(MPI_COMM_WORLD, parser);
    //	else if( parser("-sim").asString() == "sic" )
    //		sim = new Test_SICMPI(MPI_COMM_WORLD, parser);
    else if( parser("-sim").asString() == "cloud" )
        sim = new Test_CloudMPI<GridMPI_t,TFlowStep>(MPI_COMM_WORLD, parser);
    else if( parser("-sim").asString() == "cloud-acoustic" )
        sim = new Test_CloudAcousticMPI<GridMPI_t,TFlowStep>(MPI_COMM_WORLD, parser);
    else if( parser("-sim").asString() == "siccloud" )
        sim = new Test_SICCloudMPI<GridMPI_t,TFlowStep>(MPI_COMM_WORLD, parser);
    // else if( parser("-sim").asString() == "2Dsbi" )
    //     sim = new Test_2DSBIMPI(MPI_COMM_WORLD, parser);
    else if( parser("-sim").asString() == "cha" )
        sim = new Test_ChannelFlowMPI<GridMPI_t,TFlowStep>(MPI_COMM_WORLD, parser);
    // else if( parser("-sim").asString() == "hit" )
    //     sim = new Test_HITMPI(MPI_COMM_WORLD, parser);
// //    else if( parser("-sim").asString() == "hitcloud" )
// //        sim = new Test_HITCloudMPI(MPI_COMM_WORLD, parser);
    // else if( parser("-sim").asString() == "tgv" )
    //     sim = new Test_TGVMPI(MPI_COMM_WORLD, parser);
// //    else if( parser("-sim").asString() == "kh" )
// //        sim = new Test_KHinstabilityMPI(MPI_COMM_WORLD, parser);
    // else if( parser("-sim").asString() == "osc" )
    //     sim = new Test_OscillatingShapeMPI(MPI_COMM_WORLD, parser);
#ifdef _USE_FFTW_
    // else if( parser("-sim").asString() == "hit" )
    //     sim = new Test_HITMPI(MPI_COMM_WORLD, parser);
#endif
    else
		if (isroot)
		{
			printf("-sim value not recognized. Aborting.\n");
			abort();
		}
		else abort();

	sim->setup();

	double wallclock;

	{
		Timer timer;

		timer.start();

#ifdef _QPX_
		HPM_Start("Simulation");
#endif
		sim->run();
#ifdef _QPX_
		HPM_Stop("Simulation");
#endif
		wallclock = timer.stop();
	}

	// sim->dispose();

	delete sim;

	sim = NULL;

	MPI_Finalize();	// peh: dropped C++ bindings

	return 0;
}
