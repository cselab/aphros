/*
 *  Test_SteadyStateMPI.h
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_STEADYSTATEMPI_H_MT52JRKU
#define TEST_STEADYSTATEMPI_H_MT52JRKU

#include <mpi.h>

#include "Test_SteadyState.h"
#include "GridMPI.h"
#include "BlockLabMPI.h"
#include "OutputProcessingMPI.h"


template <typename TGrid, typename TStepper, template <typename> class TSlice=SliceMPI>
class Test_SteadyStateMPI: public virtual Test_SteadyState<TGrid,TStepper,TSlice>
{
public:

    Test_SteadyStateMPI(const MPI_Comm comm, ArgumentParser& P) :
        Test_SteadyState<TGrid,TStepper,TSlice>(P), m_comm_world(comm)
    {
        int rank;
        MPI_Comm_rank(m_comm_world, &rank);
        this->isroot = (0 == rank);
        if (!this->isroot) this->VERBOSITY = 0;
    }
    virtual ~Test_SteadyStateMPI() {}

    virtual void setup();


protected:
    MPI_Comm m_comm_world;

    int XPESIZE, YPESIZE, ZPESIZE;

    // case specific
    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////      TEST STEADY STATE   MPI    ///////////////\n");
        printf("////////////////////////////////////////////////////////////\n");
        typedef typename TGrid::BlockType B;
        std::cout << "Domain size:   [" << this->BPDX*B::sizeX*XPESIZE;
        std::cout << " x " << this->BPDY*B::sizeY*YPESIZE;
        std::cout << " x " <<  this->BPDZ*B::sizeZ*ZPESIZE << "]" << std::endl;

        std::cout << "Domain extent: [" << Simulation_Environment::extents[0];
        std::cout << " x " << Simulation_Environment::extents[1];
        std::cout << " x " <<  Simulation_Environment::extents[2] << "]" << std::endl;
    }
    virtual void _setup_parameter();

    virtual void _init()
    {
        Test_SteadyState<TGrid,TStepper,TSlice>::_init();
        MPI_Barrier(m_comm_world);
    }
};


// class implementation
template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_setup_parameter()
{
    Test_SteadyState<TGrid,TStepper,TSlice>::_setup_parameter();

    XPESIZE = this->parser("-xpesize").asInt(2);
    YPESIZE = this->parser("-ypesize").asInt(2);
    ZPESIZE = this->parser("-zpesize").asInt(2);

    const int bpdx = this->BPDX;
    const int bpdy = this->BPDY;
    const int bpdz = this->BPDZ;
    const int BPD_PE_MAX = 
        std::max(std::max(bpdx*XPESIZE, bpdy*YPESIZE), bpdz*ZPESIZE);
    Simulation_Environment::extents[0] = 
        Simulation_Environment::extent * 
        (bpdx*XPESIZE) / static_cast<double>(BPD_PE_MAX);
    Simulation_Environment::extents[1] = 
        Simulation_Environment::extent * 
        (bpdy*YPESIZE) / static_cast<double>(BPD_PE_MAX);
    Simulation_Environment::extents[2] = 
        Simulation_Environment::extent * 
        (bpdz*ZPESIZE) / static_cast<double>(BPD_PE_MAX);
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyStateMPI<TGrid,TStepper,TSlice>::setup()
{
    _setup_parameter();

    if (this->isroot)
        _print_case_header();

    this->grid = new TGrid(
        XPESIZE, YPESIZE, ZPESIZE, 
        this->BPDX, this->BPDY, this->BPDZ, 
        Simulation_Environment::extent, m_comm_world);
    assert(this->grid != NULL);

    this->stepper = new TStepper(*(this->grid), this->parser, this->VERBOSITY);
    assert(this->stepper != NULL);

    this->dumper = new OutputProcessingMPI<TGrid,TSlice>(this->parser, *(this->grid), this->isroot);
    assert(this->dumper != NULL);
    this->dumper->register_all(*(this->grid));

    const std::string path = this->parser("-fpath").asString(".");
    this->_setup_ic();
}

#endif /* TEST_STEADYSTATEMPI_H_MT52JRKU */
