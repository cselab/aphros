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
#ifdef _QPXEMU_
#include <xmmintrin.h>
#endif

#ifdef __bgq__
#include <spi/include/l1p/sprefetch.h>
#endif

#include "Test_SteadyState.h"
#include "FlowStep_LSRK3.h"
// #include "Test_ShockBubble.h"
// /* #include "Test_SIC.h" */
#include "Test_Advection.h"
#include "Test_Cloud.h"
#include "Test_CloudAcoustic.h"
#include "Test_SICCloud.h"
// /* #include "Test_CVT.h" */
#include "Test_ChannelFlow.h"
// #include "Test_ChannelFlowMultiPhase.h"
// #include "Test_ShearFlow.h"
// #include "Test_TGV.h"
// #include "Test_SOD.h"

using namespace std;

Simulation * sim = NULL;

#ifdef _GLUT_VIZ
struct VisualSupport
{
	static void display(){}

	static void idle(void)
	{
		glClear(GL_COLOR_BUFFER_BIT);
		sim->run();
		glutSwapBuffers();
	}

	static void run(int argc, const char ** argv)
	{
		static bool bSetup = false;

		if (!bSetup)
		{
			setup(argc, argv);
			bSetup = true;
		}

		glutDisplayFunc(display);
		glutIdleFunc(idle);

		glutMainLoop();
	}

	static void setup(int argc,  const char ** argv)
	{
		glutInit(&argc, const_cast<char **>(argv));
		glutInitWindowSize(800,800);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH | GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );

		glutCreateWindow("MPCF-node");

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		glOrtho(0.0, 1.0, 0.0, 1.0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glEnable(GL_TEXTURE_2D);

		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	}
};
#endif

int main (int argc, const char ** argv)
{
#ifdef _GIT_HASH_
    const string git_hash(_GIT_HASH_);
#else
    const string git_hash("UNKNOWN");
#endif
    cout << "GIT COMMIT: " << git_hash << endl;
    cout << "Potential number of threads is " << omp_get_max_threads() << endl;

#ifdef _USE_HPM_
    MPI_Init(&argc, (char ***)&argv);
#endif

#ifdef _USE_NUMA_
	if (numa_available() < 0)
		printf("WARNING: The system does not support NUMA API!\n");
	else
		printf("NUMA API supported!\n");
#endif

    // we use one unified parser for clean overview.  use a refrence or pointer
    // to this instance if you need a parser in in your source that is included
    // in this main.cpp file
	ArgumentParser parser(argc, argv);

    // read configuration files, if available
    if (parser.exist("cubismConf"))
    {
        parser.readFile(parser("cubismConf").asString());
        const size_t NX = parser("NX").asInt();
        const size_t NY = parser("NY").asInt();
        const size_t NZ = parser("NZ").asInt();
        const size_t BX = parser("bpdx").asInt();
        const size_t BY = parser("bpdy").asInt();
        const size_t BZ = parser("bpdz").asInt();
        const bool Xgood = (NX == _BLOCKSIZE_*BX);
        const bool Ygood = (NY == _BLOCKSIZE_*BY);
        const bool Zgood = (NZ == _BLOCKSIZE_*BZ);
        const bool notGood = !(Xgood && Ygood && Zgood);
        if (notGood)
        {
            cout << "Sanity check failed... Check file: " << parser("cubismConf").asString() << endl;
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

    if( parser("-sim").asString() == "steady" )
        sim = new Test_SteadyState<Grid_t,FlowStep_LSRK3>(parser);
    // else if( parser("-sim").asString() == "sb" )
    //     sim = new Test_ShockBubble(parser);
    // /* else if( parser("-sim").asString() == "sic" ) */
    // /*     sim = new Test_SIC(parser); */
    else if( parser("-sim").asString() == "advection" )
        sim = new Test_Advection<Grid_t,FlowStep_LSRK3>(parser);
    else if( parser("-sim").asString() == "cloud" )
        sim = new Test_Cloud<Grid_t,FlowStep_LSRK3>(parser);
    else if( parser("-sim").asString() == "cloud-acoustic" )
        sim = new Test_CloudAcoustic<Grid_t,FlowStep_LSRK3>(parser);
    else if( parser("-sim").asString() == "siccloud" )
        sim = new Test_SICCloud<Grid_t,FlowStep_LSRK3>(parser);
    // else if( parser("-sim").asString() == "sod" )
    //     sim = new Test_SOD(parser);
    // else if( parser("-sim").asString() == "tgv" )
    //     sim = new Test_TGV(parser);
    else if( parser("-sim").asString() == "cha" )
        sim = new Test_ChannelFlow<Grid_t,FlowStep_LSRK3>(parser);
    // else if( parser("-sim").asString() == "champ" )
    //     sim = new Test_ChannelFlowMultiPhase(parser);
    // else if( parser("-sim").asString() == "shear" )
    //     sim = new Test_ShearFlow(parser);
    /* else if( parser("-sim").asString() == "cvt" ) */
    /*     sim = new Test_CVT(parser); */
    else
	{
		printf("Study case not defined!\n");
		abort();
	}

	sim->setup();

#if 0
	#pragma omp parallel
	{
		int hwthread = Kernel_PhysicalHWThreadID(); // returns 0-3 on each core
		if (hwthread == 0) {
			L1P_SetStreamPolicy(L1P_confirmed_or_dcbt);
			L1P_SetStreamDepth(stream_depth);
		}
	}
#else
#ifdef __bgq__
	#pragma omp parallel
	{
		int hwthread = Kernel_PhysicalHWThreadID(); // returns 0-3 on each core
		char *ptr;
		int stream_depth;
		int adaptive_flag;

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

			ptr = getenv("L1P_ADAPTIVE");
			if (ptr != NULL) {
				adaptive_flag = atoi(ptr);
				L1P_SetStreamDepth(adaptive_flag);
			}
		}
	}
#endif
#endif

	double wallclock;

	{
		Timer timer;

		timer.start();
#ifdef _MRAG_GLUT_VIZ
		VisualSupport::run(argc, argv);
#else
		sim->run();
#endif

		wallclock = timer.stop();
	}

#ifdef _USE_HPM_
	MPI_Finalize();
#endif

    return 0;
}



void set_dcbt()
{
#ifdef __bgq__
	int hwthread = Kernel_PhysicalHWThreadID(); // returns 0-3 on each core
	if (hwthread == 0) {
		L1P_SetStreamPolicy(L1P_confirmed_or_dcbt);
		L1P_SetStreamDepth(2);
	}
#endif
}


void set_dis()
{
#ifdef __bgq__
	int hwthread = Kernel_PhysicalHWThreadID(); // returns 0-3 on each core
	if (hwthread == 0) {
		L1P_SetStreamPolicy(L1P_stream_disable);
	}
#endif
}



void set_opt()
{
#ifdef __bgq__
	int hwthread = Kernel_PhysicalHWThreadID(); // returns 0-3 on each core
	if (hwthread == 0) {
		L1P_SetStreamPolicy(L1P_stream_optimistic);
	}
#endif
}
