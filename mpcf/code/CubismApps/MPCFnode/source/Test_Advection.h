/*
 *  Test_Advection.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger on 11/10/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_ADVECTION_H_8BFHY7CZ
#define TEST_ADVECTION_H_8BFHY7CZ

#include "Test_Cloud.h"


template <typename TGrid, typename TStepper, template <typename> class TSlice=Slice>
class Test_Advection : public virtual Test_Cloud<TGrid,TStepper,TSlice>
{
public:
    Test_Advection(ArgumentParser& P) :
        Test_SteadyState<TGrid,TStepper,TSlice>(P),
        Test_Cloud<TGrid,TStepper,TSlice>(P)
    { }
    virtual ~Test_Advection() { }


protected:
    typedef typename TGrid::BlockType TBlock;
    typedef typename TBlock::ElementType TElement;

    bool m_bWithShapes, m_bCollinear;
    double m_rho1, m_rho2, m_u1, m_u2, m_p1, m_p2;


    // case specific
    virtual void _post_save() { }
    virtual void _post_restart() { }

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////           TEST ADVECTION        ///////////////\n");
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
        Test_SteadyState<TGrid,TStepper,TSlice>::_setup_parameter();
        this->_parameter_advection();
    }

    virtual void _ic();
    void _parameter_advection();
    void _set_shape_advection_ic();
    void _set_collinear_advection_ic();
    void _set_normal_advection_ic();
};


///////////////////////////////////////////////////////////////////////////////
// CLASS IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////
template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_Advection<TGrid,TStepper,TSlice>::_ic()
{
    if (this->isroot)
        std::cout << "Advection Initial condition..." << std::endl;

    if (m_bWithShapes)
        this->_set_shape_advection_ic();
    else if (m_bCollinear)
        this->_set_collinear_advection_ic();
    else
        this->_set_normal_advection_ic();

    if (this->isroot)
        cout << "done." << endl;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_Advection<TGrid,TStepper,TSlice>::_parameter_advection()
{
    m_bWithShapes = this->parser("withshapes").asBool(false);
    m_bCollinear  = this->parser("collinear").asBool(false);
    m_rho1 = this->parser("rho1").asDouble(1000.0);
    m_rho2 = this->parser("rho2").asDouble(1.0);
    m_u1   = this->parser("u1").asDouble(0.0);
    m_u2   = this->parser("u2").asDouble(0.0);
    m_p1   = this->parser("p1").asDouble(1.0);
    m_p2   = this->parser("p2").asDouble(1.0);

    AdvectionBCData::boundaryElementVector.resize(2);
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_Advection<TGrid,TStepper,TSlice>::_set_shape_advection_ic()
{
    const double G1 = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double G2 = 1.0/(Simulation_Environment::GAMMA2-1.0);
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;

    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();

    Seed<shape> myseed = this->_make_shapes();

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
            double myextent[3] = { info.h, info.h, info.h } ;
            TBlock& b = *(TBlock*)info.ptrBlock;
            std::vector<shape> myshapes = myseed.get_shapes().size() ? myseed.retain_shapes(info.origin, myextent).get_shapes() : myseed.get_shapes();

            for(int iz=0; iz<TBlock::sizeZ; iz++)
                for(int iy=0; iy<TBlock::sizeY; iy++)
                    for(int ix=0; ix<TBlock::sizeX; ix++)
                    {
                        double pos[3];
                        info.pos(pos, ix, iy, iz);

                        const double inShape = eval(myshapes, pos);

                        TElement IC;

                        const double ar1 = (1-inShape)*this->m_rho1;
                        const double ar2 = inShape*this->m_rho2;
                        const double rhoMix = ar1 + ar2;
                        IC.alpha1rho1 = ar1;
                        IC.alpha2rho2 = ar2;

                        const double u = this->m_u2*inShape + this->m_u1*(1-inShape);
                        IC.ru = u*rhoMix;
                        IC.rv = 0.0;
                        IC.rw = 0.0;

                        const double a2 = inShape;
                        const double a1 = 1.0 - a2;
                        const double pressure  = this->m_p2*inShape + this->m_p1*(1-inShape);
                        const double gmixm1Inv = a1*G1 + a2*G2;
                        const double pcmix     = a1*G1*F1 + a2*G2*F2;
                        const double Ekin      = 0.5*rhoMix*u*u;
                        IC.energy = gmixm1Inv*pressure + pcmix + Ekin;
                        IC.alpha2 = a2;
                        IC.dummy = 0.0;

                        b(ix,iy,iz) = IC;
                    }
        }
    }

    // set Dirichlet inflow
    AdvectionBCData::boundaryElementVector[0].alpha1rho1 = m_rho1;
    AdvectionBCData::boundaryElementVector[0].alpha2rho2 = 0.0;
    AdvectionBCData::boundaryElementVector[0].ru = m_rho1*m_u1;
    AdvectionBCData::boundaryElementVector[0].rv = 0.0;
    AdvectionBCData::boundaryElementVector[0].rw = 0.0;
    AdvectionBCData::boundaryElementVector[0].energy = m_p1*G1 + F1*G1 + 0.5*m_rho1*m_u1*m_u1; // v = w = 0 !
    AdvectionBCData::boundaryElementVector[0].alpha2 = 0.0;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_Advection<TGrid,TStepper,TSlice>::_set_collinear_advection_ic()
{
    const double G1 = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double G2 = 1.0/(Simulation_Environment::GAMMA2-1.0);
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;

    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();

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
            TBlock& b = *(TBlock*)info.ptrBlock;

            for(int iz=0; iz<TBlock::sizeZ; iz++)
                for(int iy=0; iy<TBlock::sizeY; iy++)
                    for(int ix=0; ix<TBlock::sizeX; ix++)
                    {
                        double pos[3];
                        info.pos(pos, ix, iy, iz);

                        double rho1, rho2, u, pressure, a2;
                        if (pos[1] < 1./3.*Simulation_Environment::extents[1])
                        {
                            rho1 = m_rho1;
                            rho2 = 0.0;
                            u = m_u1;
                            pressure = m_p1;
                            a2 = 0.0;
                        }
                        else if (pos[1] < 2./3.*Simulation_Environment::extents[1])
                        {
                            rho1 = 0.0;
                            rho2 = m_rho2;
                            u = m_u2;
                            pressure = m_p2;
                            a2 = 1.0;
                        }
                        else
                        {
                            rho1 = m_rho1;
                            rho2 = 0.0;
                            u = m_u1;
                            pressure = m_p1;
                            a2 = 0.0;
                        }

                        TElement IC;

                        const double rhoMix = rho1 + rho2;
                        IC.alpha1rho1 = rho1;
                        IC.alpha2rho2 = rho2;

                        IC.ru = u*rhoMix;
                        IC.rv = 0.0;
                        IC.rw = 0.0;

                        const double a1 = 1.0 - a2;
                        const double gmixm1Inv = a1*G1 + a2*G2;
                        const double pcmix     = a1*G1*F1 + a2*G2*F2;
                        const double Ekin      = 0.5*rhoMix*u*u;
                        IC.energy = gmixm1Inv*pressure + pcmix + Ekin;
                        IC.alpha2 = a2;
                        IC.dummy = 0.0;

                        b(ix,iy,iz) = IC;
                    }
        }
    }

    // set Dirichlet inflow
    AdvectionBCData::boundaryElementVector[0].alpha1rho1 = m_rho1;
    AdvectionBCData::boundaryElementVector[0].alpha2rho2 = 0.0;
    AdvectionBCData::boundaryElementVector[0].ru = m_rho1*m_u1;
    AdvectionBCData::boundaryElementVector[0].rv = 0.0;
    AdvectionBCData::boundaryElementVector[0].rw = 0.0;
    AdvectionBCData::boundaryElementVector[0].energy = m_p1*G1 + F1*G1 + 0.5*m_rho1*m_u1*m_u1; // v = w = 0 !
    AdvectionBCData::boundaryElementVector[0].alpha2 = 0.0;

    AdvectionBCData::boundaryElementVector[1].alpha1rho1 = 0.0;
    AdvectionBCData::boundaryElementVector[1].alpha2rho2 = m_rho2;
    AdvectionBCData::boundaryElementVector[1].ru = m_rho2*m_u2;
    AdvectionBCData::boundaryElementVector[1].rv = 0.0;
    AdvectionBCData::boundaryElementVector[1].rw = 0.0;
    AdvectionBCData::boundaryElementVector[1].energy = m_p2*G2 + F2*G2 + 0.5*m_rho2*m_u2*m_u2; // v = w = 0 !
    AdvectionBCData::boundaryElementVector[1].alpha2 = 1.0;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_Advection<TGrid,TStepper,TSlice>::_set_normal_advection_ic()
{
    const double G1 = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double G2 = 1.0/(Simulation_Environment::GAMMA2-1.0);
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;

    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();

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
            TBlock& b = *(TBlock*)info.ptrBlock;

            for(int iz=0; iz<TBlock::sizeZ; iz++)
                for(int iy=0; iy<TBlock::sizeY; iy++)
                    for(int ix=0; ix<TBlock::sizeX; ix++)
                    {
                        double pos[3];
                        info.pos(pos, ix, iy, iz);

                        double rho1, rho2, u, pressure, a2;
                        if (pos[0] < 0.5*Simulation_Environment::extents[0])
                        {
                            rho1 = m_rho1;
                            rho2 = 0.0;
                            u = m_u1;
                            pressure = m_p1;
                            a2 = 0.0;
                        }
                        else
                        {
                            rho1 = 0.0;
                            rho2 = m_rho2;
                            u = m_u2;
                            pressure = m_p2;
                            a2 = 1.0;
                        }

                        TElement IC;

                        const double rhoMix = rho1 + rho2;
                        IC.alpha1rho1 = rho1;
                        IC.alpha2rho2 = rho2;

                        IC.ru = u*rhoMix;
                        IC.rv = 0.0;
                        IC.rw = 0.0;

                        const double a1 = 1.0 - a2;
                        const double gmixm1Inv = a1*G1 + a2*G2;
                        const double pcmix     = a1*G1*F1 + a2*G2*F2;
                        const double Ekin      = 0.5*rhoMix*u*u;
                        IC.energy = gmixm1Inv*pressure + pcmix + Ekin;
                        IC.alpha2 = a2;
                        IC.dummy = 0.0;

                        b(ix,iy,iz) = IC;
                    }
        }
    }

    // set Dirichlet inflow
    AdvectionBCData::boundaryElementVector[0].alpha1rho1 = m_rho1;
    AdvectionBCData::boundaryElementVector[0].alpha2rho2 = 0.0;
    AdvectionBCData::boundaryElementVector[0].ru = m_rho1*m_u1;
    AdvectionBCData::boundaryElementVector[0].rv = 0.0;
    AdvectionBCData::boundaryElementVector[0].rw = 0.0;
    AdvectionBCData::boundaryElementVector[0].energy = m_p1*G1 + F1*G1 + 0.5*m_rho1*m_u1*m_u1; // v = w = 0 !
    AdvectionBCData::boundaryElementVector[0].alpha2 = 0.0;
}

#endif /* TEST_ADVECTION_H_8BFHY7CZ */
