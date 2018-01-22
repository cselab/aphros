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

    virtual void _save()
    {
        Test_SteadyState<TGrid,TStepper,TSlice>::_save();
        MPI_Barrier(m_comm_world);
    }

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

    // MPI stuff
    XPESIZE = this->parser("-xpesize").asInt(2);
    YPESIZE = this->parser("-ypesize").asInt(2);
    ZPESIZE = this->parser("-zpesize").asInt(2);

#ifdef _NONUNIFORM_BLOCK_
    Simulation_Environment::extents[0] = this->parser("-extentX").asDouble(1.0);
    Simulation_Environment::extents[1] = this->parser("-extentY").asDouble(1.0);
    Simulation_Environment::extents[2] = this->parser("-extentZ").asDouble(1.0);
#else
    const int bpdx = this->BPDX;
    const int bpdy = this->BPDY;
    const int bpdz = this->BPDZ;
    const int BPD_PE_MAX = std::max(std::max(bpdx*XPESIZE, bpdy*YPESIZE), bpdz*ZPESIZE);
    Simulation_Environment::extents[0] = Simulation_Environment::extent*(bpdx*XPESIZE)/static_cast<double>(BPD_PE_MAX);
    Simulation_Environment::extents[1] = Simulation_Environment::extent*(bpdy*YPESIZE)/static_cast<double>(BPD_PE_MAX);
    Simulation_Environment::extents[2] = Simulation_Environment::extent*(bpdz*ZPESIZE)/static_cast<double>(BPD_PE_MAX);
#endif /* _NONUNIFORM_BLOCK_ */

    if (this->parser("restart_format").asString("zbin") == "zbin") // we prioritize zbin
        this->_serializationKernel = &DumpZBin_MPI<TGrid, StreamerSerialization>;
#ifdef _USE_HDF_
    else if (this->parser("restart_format").asString("zbin") == "h5")
        this->_serializationKernel = &DumpHDF5_MPI<TGrid, StreamerSerialization>;
#endif
    else
    {
        if (this->isroot)
            std::cerr << "ERROR: No suitable restart format chosen!" << std::endl;
        abort();
    }
    assert(this->_serializationKernel != NULL);
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyStateMPI<TGrid,TStepper,TSlice>::setup()
{
    _setup_parameter();

    if (this->isroot)
        _print_case_header();

#ifdef _NONUNIFORM_BLOCK_
    const double eX = Simulation_Environment::extents[0];
    const double eY = Simulation_Environment::extents[1];
    const double eZ = Simulation_Environment::extents[2];

    const int NBlocksX = XPESIZE * this->BPDX;
    const int NBlocksY = YPESIZE * this->BPDY;
    const int NBlocksZ = ZPESIZE * this->BPDZ;
#ifdef _WENO3_
    this->m_nonuniform = new NonUniformScheme<typename TGrid::BlockType, Weno3Coefficients>(0.0,eX,0.0,eY,0.0,eZ, NBlocksX, NBlocksY, NBlocksZ);
#else
    this->m_nonuniform = new NonUniformScheme<typename TGrid::BlockType, Weno5Coefficients_Coralic>(0.0,eX,0.0,eY,0.0,eZ, NBlocksX, NBlocksY, NBlocksZ);
#endif /* _WENO3_ */
    assert(this->m_nonuniform != NULL);

    // initialize scheme
    MeshDensityFactory mk(this->parser);
    this->m_nonuniform->init(mk.get_mesh_kernel(0), mk.get_mesh_kernel(1), mk.get_mesh_kernel(2));

    this->grid = new TGrid(&(this->m_nonuniform->get_map_x()), &(this->m_nonuniform->get_map_y()), &(this->m_nonuniform->get_map_z()),
            XPESIZE, YPESIZE, ZPESIZE, this->BPDX, this->BPDY, this->BPDZ,
            m_comm_world);
    assert(this->grid != NULL);

    vector<BlockInfo> infos = this->grid->getBlocksInfo();
    this->m_nonuniform->setup_coefficients(infos);

    this->m_nonuniform->print_mesh_statistics(this->isroot);

#else

    this->grid = new TGrid(XPESIZE, YPESIZE, ZPESIZE, this->BPDX, this->BPDY, this->BPDZ, Simulation_Environment::extent, m_comm_world);
    assert(this->grid != NULL);
#endif /* _NONUNIFORM_BLOCK_ */

    this->stepper = new TStepper(*(this->grid), this->CFL, this->parser, this->VERBOSITY);
    assert(this->stepper != NULL);

#ifdef _NONUNIFORM_BLOCK_
    this->stepper->set_hmin(this->m_nonuniform->minimum_cell_width());
#endif /* _NONUNIFORM_BLOCK_ */

    this->dumper = new OutputProcessingMPI<TGrid,TSlice>(this->parser, *(this->grid), this->isroot);
    assert(this->dumper != NULL);
    this->dumper->register_all(*(this->grid));

    const std::string path = this->parser("-fpath").asString(".");
    typename Test_SteadyState<TGrid,TStepper,TSlice>::TDeserializer _deserializer = NULL;
    if (exist_restart("zbin", path)) // we prioritize zbin (MPI topology must be the same as at serialization time)
        _deserializer = &ReadZBin_MPI<TGrid, StreamerSerialization>;
#ifdef _USE_HDF_
    else if (exist_restart("h5", path))
        _deserializer = &ReadHDF5_MPI<TGrid, StreamerSerialization>;
#endif

    this->_setup_ic(_deserializer);
}

#endif /* TEST_STEADYSTATEMPI_H_MT52JRKU */
