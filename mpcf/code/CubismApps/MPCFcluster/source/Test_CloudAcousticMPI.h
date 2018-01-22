/*
 *  Test_CloudMPI.h
 *  MPCFcluster
 *
 *  Created by Babak Hejazialhosseini on 2/25/13.
 *  Completly revised by Ursula Rasthofer in 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_CLOUDACOUSTICMPI_H_QFKFI7VT
#define TEST_CLOUDACOUSTICMPI_H_QFKFI7VT

#include "Test_CloudAcoustic.h"
#include "Test_SteadyStateMPI.h"
#include "StatisticsMPI_MultiPhase.h"

template <typename TGrid, typename TStepper, template <typename> class TSlice=SliceMPI>
class Test_CloudAcousticMPI :
    public Test_CloudAcoustic<TGrid,TStepper,TSlice>,
    public Test_SteadyStateMPI<TGrid,TStepper,TSlice>
{
public:
    Test_CloudAcousticMPI(MPI_Comm comm, ArgumentParser& _parser) :
        Test_SteadyState<TGrid,TStepper,TSlice>(_parser),
        Test_CloudAcoustic<TGrid,TStepper,TSlice>(_parser),
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>(comm, _parser)
    { }
    virtual ~Test_CloudAcousticMPI() {}


protected:
    typedef typename TGrid::BlockType TBlock;
    typedef typename TBlock::ElementType TElement;

    // case specific
    virtual void _post_save()
    {
#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
        BGC::bgc_save(this->restart_id, this->m_comm_world);
#endif
    }

    virtual void _post_restart()
    {
#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
        BGC::bgc_restart(this->restart_id, this->m_comm_world);
#endif
    }

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////     TEST CLOUD ACOUSTIC MPI     ///////////////\n");
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
        Test_CloudAcoustic<TGrid,TStepper,TSlice>::_setup_parameter();
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_setup_parameter();
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

#endif /* TEST_CLOUDACOUSTICMPI_H_QFKFI7VT */

