/*
 *  Test_SICCloudMPI.h
 *  MPCFcluster
 *
 *  Created by Fabian Wermelinger on 04/04/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_SICCLOUDMPI_H_OTVTKOMY
#define TEST_SICCLOUDMPI_H_OTVTKOMY

#include "Test_SICCloud.h"
#include "Test_SteadyStateMPI.h"
#include "StatisticsMPI_MultiPhase.h"

template <typename TGrid, typename TStepper, template <typename> class TSlice=SliceMPI>
class Test_SICCloudMPI :
    public Test_SICCloud<TGrid,TStepper,TSlice>,
    public Test_SteadyStateMPI<TGrid,TStepper,TSlice>
{
public:
    Test_SICCloudMPI(MPI_Comm comm, ArgumentParser& _parser) :
        Test_SteadyState<TGrid,TStepper,TSlice>(_parser),
        Test_SICCloud<TGrid,TStepper,TSlice>(_parser),
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>(comm, _parser)
    { }
    virtual ~Test_SICCloudMPI() {}


protected:
    // case specific
    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("//////////// TEST SHOCK INDUCED COLLAPSE MPI ///////////////\n");
        printf("////////////////////////////////////////////////////////////\n");
        typedef typename TGrid::BlockType B;
        std::cout << "Domain size:   [" << this->BPDX*B::sizeX*this->XPESIZE;
        std::cout << " x " << this->BPDY*B::sizeY*this->YPESIZE;
        std::cout << " x " <<  this->BPDZ*B::sizeZ*this->ZPESIZE << "]" << std::endl;

        std::cout << "Domain extent: [" << Simulation_Environment::extents[0];
        std::cout << " x " << Simulation_Environment::extents[1];
        std::cout << " x " <<  Simulation_Environment::extents[2] << "]" << std::endl;
    }
    virtual void _setup_parameter()
    {
        Test_SICCloud<TGrid,TStepper,TSlice>::_setup_parameter();
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_setup_parameter();
    }

    virtual void _init()
    {
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_init();
        this->_init_sic_cloud();
    }

    virtual void _analysis()
    {
        if (this->ANALYSISPERIOD != 0 && this->step_id%this->ANALYSISPERIOD == 0 && !this->bRESTART)
        {
            this->profiler.push_start("ANALYSIS");
            MultiPhaseStatistics::dumpStatistics(*this->grid, this->step_id, this->t, this->dt, this->parser);
            this->profiler.pop_stop();
        }
    }
};

#endif /* TEST_SICCLOUDMPI_H_OTVTKOMY */
