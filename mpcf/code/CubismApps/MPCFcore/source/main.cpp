/*
 *  main.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#if defined(_QPXEMU_) && _BLOCKSIZE_%4!=0
#error BLOCKSIZE NOT GOOD FOR QPX
#endif

#include <iostream>
#include <string>
#include <omp.h>

#ifdef _USE_HPM_
#include <mpi.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
//extern "C" void HPM_Init(void);
//extern "C" void HPM_Print(void);
#else
#define HPM_Start(x)
#define HPM_Stop(x)
#define MPI_Init(a,b)
#define MPI_Finalize()
#endif

#include <unistd.h>
#include <typeinfo>
#include <ArgumentParser.h>

//#include "TestTypes.h"

#include "Test_Kernel.h"
//#include "Test_LocalKernel.h"
#include "Test_LocalKernel_Val_5eq.h"
#include "Test_LocalKernel_Block_5eq.h"

//#include "Convection_CPP.h"
#include "Convection_CPP_HLLC_5eq.h"
#include "InterfaceSharpening_CPP_5eq.h"
#include "InterfaceSharpening_CPP_cell_compact_5eq.h"
//#include "Source_CPP.h"


#if defined(_QPX_) || defined(_QPXEMU_)
#include "Convection_QPX_HLLC_5eq.h"
#include "Update_QPX.h"
#include "MaxSpeedOfSound_QPX_5eq.h"
#endif


#include "Update.h"
//#include "MaxSpeedOfSound.h"
#include "MaxSpeedOfSound_CPP_5eq.h"
#include "MaxInterfaceVel_CPP_5eq.h"

using namespace std;

struct TestInfo
{
  string name;
  bool profiling, performance, accuracy;
  Real peakperf, peakbandwidth;
  int nofblocks, noftimes;

  TestInfo(string name): name(name){}
};

template<typename Test , typename Kernel> void testing(Test test, Kernel kernel, const TestInfo info)
{
  if (info.name == "all")
    printKernelName(typeid(kernel).name());
  else
    printKernelName(info.name);

  if(info.accuracy)
    test.accuracy(kernel);
  if(info.profiling)
    test.profile(kernel, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);

   printEndKernelTest();
}

template<typename Test , typename Kernel1, typename Kernel2> void comparing(Test test, Kernel1 kernel1, Kernel2 kernel2, const TestInfo info)
{
  if (info.name == "all")
    printKernelName(typeid(kernel1).name());
  else
    printKernelName(info.name);

  if (info.accuracy)
    test.accuracy(kernel2, kernel1);
  if (info.performance)
    test.performance(kernel2, kernel1, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);

  printEndKernelTest();
}


int main (int argc, const char ** argv)
{

#ifdef _GIT_HASH_
    const string git_hash(_GIT_HASH_);
#else
    const string git_hash("UNKNOWN");
#endif
    cout << "GIT COMMIT: " << git_hash << endl;

#ifdef _USE_HPM_
    MPI_Init(&argc, const_cast<char ***>(&argv));
#endif

	ArgumentParser parser(argc, argv);
        if (parser.exist("simConf"))
          parser.readFile(parser("simConf").asString());
        // print all run gonfiguration settings to stdout
        cout << "----------------------------------------------------------------------" << endl;
        cout << "* RUNTIME ARGUMENTS:" << endl;
        cout << "----------------------------------------------------------------------" << endl;
        parser.print_args();
        cout << "----------------------------------------------------------------------" << endl;

	//enable/disable the handling of denormalized numbers
#ifdef _QPXEMU_
	if (parser("-f2z").asBool(false))
#pragma omp parallel
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	}
#endif

	//kernel name
	const string kernel = parser("-kernel").asString("all");

	TestInfo info(kernel);

        //result check
        info.accuracy = parser("-accuracy").asBool(true);

        info.performance = parser("-performance").asBool(true);

	//enable/disable the performance comparison betweens kernels and their baseline
	info.profiling = parser("-profile").asBool(false);

	//nominal peak performance per core
	info.peakperf = parser("-pp").asDouble(21);

	//nominal peak bandwidth per core
	info.peakbandwidth = parser("-pb").asDouble(4.5);

	//memory footprint per thread, in terms of blocks.
	info.nofblocks = parser("-nblocks").asInt(50);

	//number of times that the kernel will be executed
	info.noftimes =  parser("-n").asInt(50);

        if (kernel=="all" and info.accuracy==true)
          std::cout << "##### WARNING: you should not run all accuracy tests at once as different reference solultions, specified as parser arguments, are required!" << std::endl;

        // get input from parser
        const Real a = parser("-lsrk-a").asDouble(0.25);
        const Real b = parser("-lsrk-b").asDouble(0.25);
        const Real dtinvh = parser("-lsrk-dtinvh").asDouble(1.5);
        const Real g1 = parser("-g1").asDouble(6.59);
        const Real g2 = parser("-g2").asDouble(1.4);
        const Real pc1 = parser("-pc1").asDouble(4.069e3);
        const Real pc2 = parser("-pc2").asDouble(0.01);
        const Real eps = parser("-eps-sharp").asDouble(0.75);
        const Real epsu0 = parser("-eps-u").asDouble(0.55);

	// C++ kernels
	// comparison to reference solution
	{
          if (kernel == "Convection_CPP_HLLC_5eq" || kernel == "all")
            testing(Test_Kernel<3>(parser), Convection_CPP_HLLC_5eq(a, dtinvh, g1, g2, pc1, pc2), info);

          if (kernel == "Diffusion_CPP_5eq" || kernel == "all")
            cout << "### WARNING: no adaption of cpp version of diffusion kernel to 5eq system, yet ###" << endl;

          if (kernel == "SurfaceTension_CPP_5eq" || kernel == "all")
            cout << "### WARNING: no adaption of cpp version of surface tension kernel to 5eq system, yet ###" << endl;

          if (kernel == "InterfaceSharpening_CPP_5eq" || kernel == "all")
            testing(Test_Kernel<1>(parser), InterfaceSharpening_CPP_5eq(a, dtinvh, eps, epsu0, g1, g2, pc1, pc2), info);

          if (kernel == "InterfaceSharpening_CPP_cell_compact_5eq" || kernel == "all")
            testing(Test_Kernel<1>(parser), InterfaceSharpening_CPP_cell_compact_5eq(a, dtinvh, eps, epsu0, g1, g2, pc1, pc2), info);


          if (kernel == "Update_CPP" || kernel == "all")
            testing(Test_LocalKernel_Block_5eq(parser), Update_CPP(b), info);

          if (kernel == "MaxSpeedOfSound_CPP_5eq" || kernel == "all")
            testing(Test_LocalKernel_Val_5eq(parser), MaxSpeedOfSound_CPP_5eq(g1, g2, pc1, pc2), info);

          if (kernel == "MaxInterfaceVel_CPP_5eq" || kernel == "all")
            testing(Test_LocalKernel_Val_5eq(parser), MaxInterfaceVel_CPP_5eq(), info);
        }

    // Vectorized
#if defined(_QPX_) || defined(_QPXEMU_)
    {
        if (kernel == "Convection_QPX_HLLC_5eq" || kernel == "all")
            testing(Test_Kernel<3>(parser), Convection_QPX_HLLC_5eq(a, dtinvh, g1, g2, pc1, pc2), info);

        if (kernel == "MaxSpeedOfSound_QPX_5eq" || kernel == "all")
            testing(Test_LocalKernel_Val_5eq(parser), MaxSpeedOfSound_QPX_5eq(g1, g2, pc1, pc2), info);

        if (kernel == "comparison_MaxSpeedOfSound_QPX_5eq" || kernel == "all")
        {
            comparing(Test_LocalKernel_Val_5eq(parser), MaxSpeedOfSound_CPP_5eq(g1, g2, pc1, pc2), MaxSpeedOfSound_QPX_5eq(g1, g2, pc1, pc2), info);
            // testing(Test_LocalKernel_Val_5eq(parser), MaxSpeedOfSound_QPX_5eq(g1, g2, pc1, pc2), info);
        }

        if (kernel == "comparison_Convection_QPX_HLLC_5eq" || kernel == "all")
            comparing(Test_Kernel<3>(parser), Convection_CPP_HLLC_5eq(a, dtinvh, g1, g2, pc1, pc2), Convection_QPX_HLLC_5eq(a, dtinvh, g1, g2, pc1, pc2), info);

        if (kernel == "mytest_L")
            comparing(Test_LocalKernel_Val_5eq(parser), MaxSpeedOfSound_CPP_5eq(g1, g2, pc1, pc2), MaxSpeedOfSound_QPX_5eq(g1, g2, pc1, pc2), info);
    }
#endif /* QPX */

// TODO: the following is only termporaty to check functionality
        {
          if (kernel == "comparison_Convection_CPP_HLLC_5eq" || kernel == "all")
            comparing(Test_Kernel<3>(parser), Convection_CPP_HLLC_5eq(a, dtinvh, g1, g2, pc1, pc2), Convection_CPP_HLLC_5eq(a, dtinvh, g1, g2, pc1, pc2), info);

          if (kernel == "comparison_Update_CPP" || kernel == "all")
            comparing(Test_LocalKernel_Block_5eq(parser), Update_CPP(3.57), Update_CPP(b), info); //this one should fail

          if (kernel == "comparison_MaxSpeedOfSound_CPP_5eq" || kernel == "all")
            comparing(Test_LocalKernel_Val_5eq(parser), MaxSpeedOfSound_CPP_5eq(g1, g2, pc1, pc2), MaxSpeedOfSound_CPP_5eq(g1, g2, pc1, pc2), info);
        }




/*
                {
                        Test_LocalKernel lt;

                        if (kernel == "Update_CPP" || kernel == "all")
                        {
                                Update_CPP update_kernel;
                                HPM_Start("Update_CPP");
                                lt.profile_update(update_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
                                HPM_Stop("Update_CPP");
                        }

                        if (kernel == "MaxSOS_CPP" || kernel == "all")
                        {
                                MaxSpeedOfSound_CPP maxsos_kernel;
                                HPM_Start("MaxSOS_CPP");
                                lt.profile_maxsos(maxsos_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
                                HPM_Stop("MaxSOS_CPP");
                        }

                        if (kernel == "MaxSOS_CPP_5eq" || kernel == "all")
                        {
                                MaxSpeedOfSound_CPP_5eq maxsos_kernel;
                                HPM_Start("MaxSOS_CPP_5eq");
                                lt.profile_maxsos(maxsos_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
                                HPM_Stop("MaxSOS_CPP_5eq");
                        }
                }
*/

        //QPX kernels
#if 0 //defined(_QPX_) || defined(_QPXEMU_)
        {
                if (kernel == "Convection_QPX" || kernel == "all")
                        testing(Test_Convection(), Convection_QPX(0, 1), info);

                Test_LocalKernel lt;

                /* if (kernel == "MaxSOS_QPX" || kernel == "all") */
                /* { */
                /*      MaxSpeedOfSound_QPX maxsos_kernel; */
                /*      MaxSpeedOfSound_CPP refkernel; */
                /*      lt.accuracy(maxsos_kernel, refkernel, info.accuracythreshold); */
                /*      HPM_Start("MaxSOS_QPX"); */
                /*      lt.profile_maxsos(maxsos_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes); */
                /*      HPM_Stop("MaxSOS_QPX"); */
                /* } */

                if (kernel == "Update_QPX" || kernel == "all")
                {
                        Update_CPP refkernel;
                        Update_QPX update_kernel;
                        lt.accuracy(update_kernel, refkernel, info.accuracythreshold);

                        HPM_Start("Update_QPX");
                        lt.profile_update(update_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
                        HPM_Stop("Update_QPX");
                }
        }
#endif

#ifdef _USE_HPM_
	MPI_Finalize();
#endif
	return 0;
}
