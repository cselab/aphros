/*
 *  Test_Cloud.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 2/24/13.
 *  Completly revised by Ursula Rasthofer in 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_CLOUDACOUSTIC_H_UDRSNGY5
#define TEST_CLOUDACOUSTIC_H_UDRSNGY5

#include <cstdlib>
#include "Test_Cloud.h"

template <typename TGrid, typename TStepper, template <typename> class TSlice=Slice>
class Test_CloudAcoustic: public Test_Cloud<TGrid,TStepper,TSlice>
{
public:
    Test_CloudAcoustic(ArgumentParser& P) :
        Test_SteadyState<TGrid,TStepper,TSlice>(P),
        Test_Cloud<TGrid,TStepper,TSlice>(P)
    { }
    virtual ~Test_CloudAcoustic() { }


protected:

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////       TEST CLOUD ACOUSTIC       ///////////////\n");
        printf("////////////////////////////////////////////////////////////\n");
        typedef typename TGrid::BlockType B;
        std::cout << "Domain size:   [" << this->BPDX*B::sizeX;
        std::cout << " x " << this->BPDY*B::sizeY;
        std::cout << " x " <<  this->BPDZ*B::sizeZ << "]" << std::endl;

        std::cout << "Domain extent: [" << Simulation_Environment::extents[0];
        std::cout << " x " << Simulation_Environment::extents[1];
        std::cout << " x " <<  Simulation_Environment::extents[2] << "]" << std::endl;
    }

    virtual void _setup_parameter()
    {
        Test_Cloud<TGrid,TStepper,TSlice>::_setup_parameter();
        this->_parameter_cloud_acoustics();
    }

    virtual void _ic()
    {
        if (this->parser.check("-single_phase"))
            this->_ic_single_phase();
        else
            this->_ic_multi_phase();
    }
    void _ic_single_phase();
    void _ic_multi_phase();

    void _parameter_cloud_acoustics();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////
template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_CloudAcoustic<TGrid,TStepper,TSlice>::_ic_single_phase()
{
    if (this->isroot)
        std::cout << "Cloud Initial condition..." << std::endl;

    const double G1 = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double G2 = 1.0/(Simulation_Environment::GAMMA2-1.0);
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;

    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();
#ifdef _NONUNIFORM_BLOCK_
    const double h = this->m_nonuniform->minimum_cell_width();
#else
    const double h = vInfo[0].h_gridpoint;
#endif /* _NONUNIFORM_BLOCK_ */


    typedef typename TGrid::BlockType TBlock;
    typedef typename TGrid::BlockType::ElementType TElement;

#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
        BGC::BPDX = this->BPDX;
        BGC::BPDY = this->BPDY;
        BGC::BPDZ = this->BPDZ;

        const int xpesize = this->parser("-xpesize").asInt(1);
        const int ypesize = this->parser("-ypesize").asInt(1);
        const int zpesize = this->parser("-zpesize").asInt(1);

        const int NBX = this->BPDX * xpesize;
        const int NBY = this->BPDY * ypesize;
        const int NBZ = this->BPDZ * zpesize;

        // loop all block and check whether they are close to the boundary
        for(int i=0; i<(int)vInfo.size(); i++)
        {
            BlockInfo info = vInfo[i];

            BGC::BoundaryBlock dummy;
            if (info.index[0] == 0)
                BGC::bgblock_list_dir0_side0.push_back(dummy);
            if (info.index[0] == (NBX-1))
                BGC::bgblock_list_dir0_side1.push_back(dummy);
            if (info.index[1] == 0)
                BGC::bgblock_list_dir1_side0.push_back(dummy);
            if (info.index[1] == (NBY-1))
                BGC::bgblock_list_dir1_side1.push_back(dummy);
            if (info.index[2] == 0)
                BGC::bgblock_list_dir2_side0.push_back(dummy);
            if (info.index[2] == (NBZ-1))
                BGC::bgblock_list_dir2_side1.push_back(dummy);
        }
#endif


#pragma omp parallel for
    for(int i=0; i<(int)vInfo.size(); i++)
    {
        BlockInfo info = vInfo[i];
        TBlock& b = *(TBlock*)info.ptrBlock;

        int s[3] = {0,0,0};
        int e[3] = {TBlock::sizeX,TBlock::sizeY,TBlock::sizeZ};

#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
        int delta_ix = 0;
        int delta_iy = 0;
        int delta_iz = 0;
        int bnx = 3;
        int bny = 3;

        if (info.index[0] == 0)
            s[0] = -3;
        if (info.index[0] == (NBX-1))
            e[0] = TBlock::sizeX+3;
        if (info.index[1] == 0)
            s[1] = -3;
        if (info.index[1] == (NBY-1))
            e[1] = TBlock::sizeY+3;
        if (info.index[2] == 0)
            s[2] = -3;
        if (info.index[2] == (NBZ-1))
            e[2] = TBlock::sizeZ+3;
#endif

        for(int iz=s[2]; iz<e[2]; iz++)
            for(int iy=s[1]; iy<e[1]; iy++)
                for(int ix=s[0]; ix<e[0]; ix++)
                {
                    TElement myvalues;
                    myvalues.alpha1rho1= CloudData::rho1;
                    myvalues.alpha2rho2= 0.0;
                    myvalues.ru = 0.0;
                    myvalues.rv = 0.0;
                    myvalues.rw = 0.0;
                    myvalues.energy = G1*CloudData::p1 + F1*G1;
                    myvalues.alpha2 = 0.0;

#if !defined(_CHARACTERISTIC_1D_BOUNDARY_)
                    b(ix,iy,iz) = myvalues;
#else
                    if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX) // regular cell in block
                        b(ix,iy,iz) = myvalues;
                    else
                        BGC::bgc_set_block<TBlock>(ix,iy,iz,info,myvalues);
#endif
                }
    }

    if (this->isroot)
        std::cout << "done." << endl;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_CloudAcoustic<TGrid,TStepper,TSlice>::_ic_multi_phase()
{
    if (this->isroot)
        std::cout << "Cloud Initial condition..." << std::endl;

    const double G1 = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double G2 = 1.0/(Simulation_Environment::GAMMA2-1.0);
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;

    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();
#ifdef _NONUNIFORM_BLOCK_
    const double h = this->m_nonuniform->minimum_cell_width();
#else
    const double h = vInfo[0].h_gridpoint;
#endif /* _NONUNIFORM_BLOCK_ */

    Seed<shape> myseed = this->_make_shapes();

    typedef typename TGrid::BlockType TBlock;
    typedef typename TGrid::BlockType::ElementType TElement;

#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
        BGC::BPDX = this->BPDX;
        BGC::BPDY = this->BPDY;
        BGC::BPDZ = this->BPDZ;

        const int xpesize = this->parser("-xpesize").asInt(1);
        const int ypesize = this->parser("-ypesize").asInt(1);
        const int zpesize = this->parser("-zpesize").asInt(1);

        const int NBX = this->BPDX * xpesize;
        const int NBY = this->BPDY * ypesize;
        const int NBZ = this->BPDZ * zpesize;

        // loop all block and check whether they are close to the boundary
        for(int i=0; i<(int)vInfo.size(); i++)
        {
            BlockInfo info = vInfo[i];

            BGC::BoundaryBlock dummy;
            if (info.index[0] == 0)
                BGC::bgblock_list_dir0_side0.push_back(dummy);
            if (info.index[0] == (NBX-1))
                BGC::bgblock_list_dir0_side1.push_back(dummy);
            if (info.index[1] == 0)
                BGC::bgblock_list_dir1_side0.push_back(dummy);
            if (info.index[1] == (NBY-1))
                BGC::bgblock_list_dir1_side1.push_back(dummy);
            if (info.index[2] == 0)
                BGC::bgblock_list_dir2_side0.push_back(dummy);
            if (info.index[2] == (NBZ-1))
                BGC::bgblock_list_dir2_side1.push_back(dummy);
        }
#endif


#pragma omp parallel
    {
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif

#pragma omp for
        for(int i=0; i<(int)vInfo.size(); i++)
        {
            BlockInfo info = vInfo[i];
            double mystart[3] = {info.origin[0], info.origin[1], info.origin[2]};
            double myextent[3] = { info.block_extent[0], info.block_extent[1], info.block_extent[2] } ;

            TBlock& b = *(TBlock*)info.ptrBlock;

            int s[3] = {0,0,0};
            int e[3] = {TBlock::sizeX,TBlock::sizeY,TBlock::sizeZ};

#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
            int delta_ix = 0;
            int delta_iy = 0;
            int delta_iz = 0;
            int bnx = 3;
            int bny = 3;

            if (info.index[0] == 0)
                s[0] = -3;
            if (info.index[0] == (NBX-1))
                e[0] = TBlock::sizeX+3;
            if (info.index[1] == 0)
                s[1] = -3;
            if (info.index[1] == (NBY-1))
                e[1] = TBlock::sizeY+3;
            if (info.index[2] == 0)
                s[2] = -3;
            if (info.index[2] == (NBZ-1))
                e[2] = TBlock::sizeZ+3;
#endif

            std::vector<shape> myshapes = myseed.get_shapes();

            for(int iz=s[2]; iz<e[2]; iz++)
                for(int iy=s[1]; iy<e[1]; iy++)
                    for(int ix=s[0]; ix<e[0]; ix++)
                    {
                        double p[3];
                        info.pos(p, ix, iy, iz);

                        // zero initial velocity is assumed
                        TElement myvalues = get_ic_conserved<TElement>(p, h, myshapes, this->m_clData);

#if !defined(_CHARACTERISTIC_1D_BOUNDARY_)
                        b(ix,iy,iz) = myvalues;
#else
                        if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX) // regular cell in block
                            b(ix,iy,iz) = myvalues;
                        else
                            BGC::bgc_set_block<TBlock>(ix,iy,iz,info,myvalues);
#endif
                    }
        }
    }

    if (this->m_clData.laplacep)
    {
        if (this->isroot)
            std::cout << "pressure relaxation for initial field currenrly not supported: ask FABIAN for further details" << std::endl;
        // _relax();
        // _set_energy(*grid);
    }

    if (this->isroot)
        std::cout << "done." << endl;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_CloudAcoustic<TGrid,TStepper,TSlice>::_parameter_cloud_acoustics()
{
    Lab dummy;
    if (dummy.name() != std::string("BlockLabCloud_1DCharNonReflect_AcousticForcing") &&
            dummy.name() != std::string("BlockLabCloud_Dirichlet_AcousticForcing"))
    {
        if (this->isroot)
            std::cout << "ERROR: Test_CloudAcoustic: Wrong BlockLab!" << std::endl;
        std::abort();
    }

#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
    BGC::pamb = CloudData::p1;
    BGC::L = this->parser("-L-boundary").asDouble(1.0);
    BGC::lambda = this->parser("-lambda").asDouble(0.75);
#endif

    // this data is defined in TestLabs.h
    CloudDataAcoustic::p_amplitude = this->parser("p_amplitude").asDouble(10.0); // bar
    CloudDataAcoustic::frequency   = this->parser("frequency").asDouble(295000.0); // Hz
    CloudDataAcoustic::t0   = this->parser("signal_t0").asDouble(0.);  // [time]
    CloudDataAcoustic::phase0   = this->parser("signal_phase0").asDouble(0.); // [-]
    CloudDataAcoustic::sigma   = this->parser("signal_sigma").asDouble(0.);  // [-]

    CloudDataAcoustic::boundaryElement0.alpha1rho1 = CloudData::rho1;
    CloudDataAcoustic::boundaryElement0.alpha2rho2 = 0.0;
    CloudDataAcoustic::boundaryElement0.ru         = 0.0;
    CloudDataAcoustic::boundaryElement0.rv         = 0.0;
    CloudDataAcoustic::boundaryElement0.rw         = 0.0;
    CloudDataAcoustic::boundaryElement0.energy = 0.0; // is set in the lab
    CloudDataAcoustic::boundaryElement0.alpha2 = 0.0;
}

#endif /* TEST_CLOUDACOUSTIC_H_UDRSNGY5 */
