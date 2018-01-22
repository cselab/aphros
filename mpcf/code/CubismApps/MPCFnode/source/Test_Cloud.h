/*
 *  Test_Cloud.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 2/24/13.
 *  Completly revised by Ursula Rasthofer in 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_CLOUD_H_XEYR5OG6
#define TEST_CLOUD_H_XEYR5OG6

#include <cassert>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "Test_SteadyState.h"
#include "Shapes.h"

struct CloudICData
{
    double R0;
    double extent[3];
    double epsilon;
    double prox;
    double buff;
    bool laplacep;
    bool integral;
};


template <typename TGrid, typename TStepper, template <typename> class TSlice=Slice>
class Test_Cloud: public virtual Test_SteadyState<TGrid,TStepper,TSlice>
{
public:
    Test_Cloud(ArgumentParser& P) :
        Test_SteadyState<TGrid,TStepper,TSlice>(P)
    { }
    virtual ~Test_Cloud() { }


protected:
    typedef typename TGrid::BlockType TBlock;
    typedef typename TBlock::ElementType TElement;

    CloudICData m_clData;

    // case specific
    virtual void _post_save()
    {
        if (this->isroot)
            std::cout << "ODE BC serialization currently not supported on node layer" << std::endl;
    }
    virtual void _post_restart()
    {
        if (this->isroot)
            std::cout << "ODE BC serialization currently not supported on node layer" << std::endl;
    }

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////             TEST CLOUD          ///////////////\n");
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
        this->_parameter_cloud();
    }

    virtual void _ic();

    virtual Seed<shape> _make_shapes()
    {
        std::string b_file = this->parser("-cloud-dat").asString("cloud.dat");
        Seed<shape> myseed;
        myseed.make_shapes(b_file, this->grid->getH(), this->isroot);
        return myseed;
    }

    virtual void _relax() { if (this->isroot) std::cout << "Pressure relaxation is not implemented on the node layer" << std::endl; }
    virtual void _set_energy()
    {
        std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();

        const double g1 = Simulation_Environment::GAMMA1;
        const double g2 = Simulation_Environment::GAMMA2;
        const double pc1 = Simulation_Environment::PC1;
        const double pc2 = Simulation_Environment::PC2;
        const double g1m1Inv = 1.0/(g1-1.0);
        const double g2m1Inv = 1.0/(g2-1.0);

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
                            //(fabianw; Fri 23 Oct 2015 10:44:36 AM CEST) WARNING: assumes ke = 0
                            // energy field contains pressure as a result of pressure relaxation

                            const double alpha2 = b(ix,iy,iz).alpha2;
                            const double alpha1 = 1.0-alpha2;
                            const double gmix_m1Inv = alpha1*g1m1Inv + alpha2*g2m1Inv;
                            const double pcmix = alpha1*g1*pc1*g1m1Inv + alpha2*g2*pc2*g2m1Inv;
                            b(ix,iy,iz).energy = gmix_m1Inv*b(ix,iy,iz).energy + pcmix;
                        }
            }
        }
    }

    void _parameter_cloud();
};

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////
// small helpers
static inline double _get_energy(const double pressure, const double alpha2)
{
    // REMARK: Simulation_Environment::GAMMA1, ... are stored as
    // Real, but if the are not very fancy values, this is ok
    const double g1 = Simulation_Environment::GAMMA1;
    const double g2 = Simulation_Environment::GAMMA2;
    const double pc1 = Simulation_Environment::PC1;
    const double pc2 = Simulation_Environment::PC2;
    const double g1m1Inv = 1.0/(g1-1.0);
    const double g2m1Inv = 1.0/(g2-1.0);

    const double alpha1 = 1.0-alpha2;
    const double gmix_m1Inv = alpha1*g1m1Inv + alpha2*g2m1Inv;
    const double pcmix = alpha1*g1*pc1*g1m1Inv + alpha2*g2*pc2*g2m1Inv;
    const double energy = gmix_m1Inv*pressure + pcmix;

    return energy;
}

// Variations of alpha2 and pressure initialization functions
///////////////////////////////////////////////////////////////////////////////
static inline double _alpha_std(const double p[3], const double h, const double alpha2, const CloudICData& data)
{
    return alpha2;
}

static inline double _alpha_tiwari(const double p[3], const double h, const double alpha2, const CloudICData& data)
{
#if defined(_BCLABCLOUDSYMABSORB_) || defined(_BCLABCLOUDSYM1DCHARNREF_)
    const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
#else
    const double r = std::sqrt((p[0] - 0.5*data.extent[0])*(p[0] - 0.5*data.extent[0]) + (p[1] - 0.5*data.extent[1])*(p[1] - 0.5*data.extent[1]) + (p[2] - 0.5*data.extent[2])*(p[2] - 0.5*data.extent[2]));
#endif
    return 0.5*(1.0 - std::tanh(0.5 * (r - data.R0)/(data.epsilon*h)));
}

static inline double _pressure_tiwari(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    // initial conditions according to Tiwari et al. 2013
    // careful: this option overwrites all fields
#if defined(_BCLABCLOUDSYMABSORB_) || defined(_BCLABCLOUDSYM1DCHARNREF_)
    const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
#else
    const double r = std::sqrt((p[0] - 0.5*data.extent[0])*(p[0] - 0.5*data.extent[0]) + (p[1] - 0.5*data.extent[1])*(p[1] - 0.5*data.extent[1]) + (p[2] - 0.5*data.extent[2])*(p[2] - 0.5*data.extent[2]));
#endif

    if (r <= data.R0)
        return CloudData::p2;
    else
        return CloudData::p1 - (data.R0/r) * (CloudData::p1 -  CloudData::p2);
}

static inline double _pressure_jump(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    // should not be used: pressure field with jump at the interface
    // causes strong oscillations in the pressure field
    return CloudData::p1*(1.0-alpha2) + CloudData::p2*alpha2;
}

template <int _COS>
static inline double _pressure_eqCloud(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    // entire cloud at constant pressure
    // increase of pressure towards the domain boundary
#if defined(_BCLABCLOUDSYMABSORB_) || defined(_BCLABCLOUDSYM1DCHARNREF_)
    const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
#else
    const double r = std::sqrt((p[0] - 0.5*data.extent[0])*(p[0] - 0.5*data.extent[0]) + (p[1] - 0.5*data.extent[1])*(p[1] - 0.5*data.extent[1]) + (p[2] - 0.5*data.extent[2])*(p[2] - 0.5*data.extent[2]));
#endif
    const double dist = r - data.R0;

    if (dist <= data.buff)
        return CloudData::p2;
    else
    {
        if (_COS)
        {
            const double alpha = M_PI*std::min(1., std::max(0., (data.prox - (dist-data.buff))/data.prox));
            return CloudData::p2 + (0.5 + 0.5 * std::cos(alpha)) * (CloudData::p1 -  CloudData::p2);
        }
        else
            return CloudData::p2 + std::tanh(2.0*dist/data.prox) * (CloudData::p1 -  CloudData::p2);
    }
}

template <int _COS>
static inline double _pressure_tiwariCloud(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    // following ideas presented in Tiwari et al 2015
    const double dist = distance(v_shapes, p);

    if (dist <= data.buff)
        return CloudData::p2;
    else
    {
        if(_COS)
        {
            const double alpha = M_PI*std::min(1., std::max(0., (data.prox - (dist-data.buff))/data.prox));
            return CloudData::p2 + (0.5 + 0.5 * std::cos(alpha)) * (CloudData::p1 -  CloudData::p2);
        }
        else
            return CloudData::p2 + std::tanh(2.0*dist/data.prox) * (CloudData::p1 -  CloudData::p2);
    }
}

// explicit instantiation
template double _pressure_eqCloud<0>(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data); // tanh
template double _pressure_eqCloud<1>(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data); // cos
template double _pressure_tiwariCloud<0>(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data); // tanh
template double _pressure_tiwariCloud<1>(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data); // cos
///////////////////////////////////////////////////////////////////////////////

// global function handles for pressure and alpha
double (*pressure_handle)(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data);
double (*alpha_handle)(const double p[3], const double h, const double alpha2, const CloudICData& data);
///////////////////////////////////////////////////////////////////////////////


// ic for cloud: primitive variables
// pressure goes to energy and velocity to momentum
// (p is position in space)
template<typename T>
inline T get_ic_primitives(const double p[3], const double h, const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    T out;
    out.clear();

    out.alpha2 = alpha_handle(p, h, alpha2, data);

    const double pres  = pressure_handle(p, h, v_shapes, alpha2, data);
    out.energy     = pres;
#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
    // store the pressure in the dummy as it is needed for characteristic-based boundary conditions
    // HINT: if dummy is need for something else we have to recompute the pressure from energy in _ic()
    out.dummy      = pres;
#endif

    return out;
}

template<typename T>
inline T get_ic_conserved(const double p[3], double h, const std::vector< shape >& v_shapes, const CloudICData& data)
{
    T out;
    out.clear();

    if (data.integral)
    {
      const double myh = h * 0.5;

      T samples[3][3][3];
      T zintegrals[3][3];
      T yzintegrals[3];

      const double x0[3] = {p[0] -myh, p[1] - myh, p[2] - myh};

      for(int iz=0; iz<3; iz++)
          for(int iy=0; iy<3; iy++)
              for(int ix=0; ix<3; ix++)
              {
                  double mypos[3] = {x0[0]+ix*myh, x0[1]+iy*myh, x0[2]+iz*myh};
                  for (int d = 0; d < 3; d++)
                    if (Simulation_Environment::BC_PERIODIC[d]) {
                      if (mypos [d] < 0.0) mypos [d] += static_cast<double>(Simulation_Environment::extents [d]);
                      if (mypos [d] > Simulation_Environment::extents[d]) mypos [d] -= static_cast<double>(Simulation_Environment::extents [d]);
                  }
                  if (v_shapes.size()>0)
                    samples[iz][iy][ix] = get_ic_primitives<T>(mypos, h, v_shapes, eval(v_shapes, mypos), data);
                  else
                    samples[iz][iy][ix] = get_ic_primitives<T>(mypos, h, v_shapes, 0, data);
              }

      for(int iy=0; iy<3; iy++)
          for(int ix=0; ix<3; ix++)
              zintegrals[iy][ix] = (1/6.) * samples[0][iy][ix]+(2./3) * samples[1][iy][ix]+(1/6.)* samples[2][iy][ix];

      for(int ix=0; ix<3; ix++)
          yzintegrals[ix] = (1./6)*zintegrals[0][ix] + (2./3)*zintegrals[1][ix]+(1./6)*zintegrals[2][ix];

      out = (1./6) * yzintegrals[0]+(2./3) * yzintegrals[1]+(1./6)* yzintegrals[2];
    }
    else
    {
      if (v_shapes.size()>0)
        out = get_ic_primitives<T>(p, h, v_shapes, eval(v_shapes, p), data);
      else
        out = get_ic_primitives<T>(p, h, v_shapes, 0, data);
    }

    // now we compute the conserved quantities
    out.alpha1rho1 = (1.0-out.alpha2)*CloudData::rho1;
    out.alpha2rho2 = out.alpha2*CloudData::rho2;
    // we assume zero velocities as also done below for the energy
//    const double rho_mix = out.alpha1rho1 + out.alpha2rho2;
//    out.ru = rho_mix*out.ru;
//    out.rv = rho_mix*out.rv;
//    out.rw = rho_mix*out.rw;
    // check whether pressure relaxtion should be performed
    // if this is the case, we keep the pressure and solve
    // for laplace p = 0; the energy is set afterwards by a
    // separate function
    if (not data.laplacep)
        out.energy = _get_energy(out.energy, out.alpha2);

   return out;
}

///////////////////////////////////////////////////////////////////////////////
// CLASS IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////
template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_Cloud<TGrid,TStepper,TSlice>::_ic()
{
    if (this->isroot)
        std::cout << "Cloud Initial condition..." << std::endl;

    const double G1 = Simulation_Environment::GAMMA1-1;
    const double G2 = Simulation_Environment::GAMMA2-1;
    const double PC1 = Simulation_Environment::PC1;
    const double PC2 = Simulation_Environment::PC2;
    const double F1 = G1*PC1;
    const double F2 = G2*PC2;

    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();
#ifdef _NONUNIFORM_BLOCK_
    const double h = this->m_nonuniform->minimum_cell_width();
#else
    const double h = vInfo[0].h_gridpoint;
#endif /* _NONUNIFORM_BLOCK_ */

    Seed<shape> myseed = this->_make_shapes();

    // rasthofer May 2016
    // set-up boundary ghost cells
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

    BGC::pamb = CloudData::p1;
    BGC::L = this->parser("-L-boundary").asDouble(1.0);
    BGC::lambda = this->parser("-lambda").asDouble(0.75);

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

    // some debug output
    /*
       {
       cout << "\n-----------------------------------------------------------\n";
       cout << "Boundary blocks:\n" ;
       cout << "xdir lhs: " << BGC::bgblock_list_dir0_side0.size() << " blocks\n";
       cout << "xdir rhs: " << BGC::bgblock_list_dir0_side1.size() << " blocks\n";
       cout << "ydir lhs: " << BGC::bgblock_list_dir1_side0.size() << " blocks\n";
       cout << "ydir rhs: " << BGC::bgblock_list_dir1_side1.size() << " blocks\n";
       cout << "zdir lhs: " << BGC::bgblock_list_dir2_side0.size() << " blocks\n";
       cout << "zdir rhs: " << BGC::bgblock_list_dir2_side1.size() << " blocks\n";
       cout << "-----------------------------------------------------------" << endl;
       }
       */
#endif

    // setting maxlayer-block=0 is fine for all ic, except for tiwari-cloud
    // in this case, each rank as to see its closest bubbles
    // recommended value for this case is 4, but you should carefully check
    // your ic before subitting a production run
    const int maxlayer = this->parser("maxlayer-block").asInt(0);

#pragma omp parallel
    {
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif

        {
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

#ifdef _KEEPALL_
                std::vector<shape> myshapes = myseed.get_shapes();
#else
                // April 20, 2016: rasthofer removed due to extended block bounding box below
                //vector<shape> myshapes = myseed.get_shapes().size() ? myseed.retain_shapes(info.origin, myextent).get_shapes() : myseed.get_shapes();
                //printf("processing %d out of %d. for this block i have: %d bubbles\n", i, (int)vInfo.size(), myshapes.size());

                std::vector<shape> myshapes;

                if (myseed.get_shapes().size() != 0)
                {
                    int bbox =1;
                    int layerplus = 0;
                    while(true)
                    {
                        double ms[3], me[3];
                        for (size_t rr=0; rr<3; rr++)
                        {
                            ms[rr] = mystart[rr] - bbox*info.block_extent[rr];
                            me[rr] = (1+2*bbox)*myextent[rr];
                        }

                        myshapes = myseed.retain_shapes(ms,me).get_shapes();

                        if (myshapes.size() == 0)
                            bbox++;
                        else
                        {
                            if (layerplus<maxlayer)
                            {
                                bbox++;
                                layerplus++;
                            }
                            else
                                break;
                        }
                    }
                }
                else
                    myshapes = myseed.get_shapes();
#endif

                {
                    for(int iz=s[2]; iz<e[2]; iz++)
                        for(int iy=s[1]; iy<e[1]; iy++)
                            for(int ix=s[0]; ix<e[0]; ix++)
                            {
                                double p[3];
                                info.pos(p, ix, iy, iz);

                                TElement myvalues;

                                myvalues  = get_ic_conserved<TElement>(p, h, myshapes, m_clData);
#if !defined(_CHARACTERISTIC_1D_BOUNDARY_)
                                b(ix,iy,iz) = myvalues;
#else
                                if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX) // regular cell in block
                                    b(ix,iy,iz) = myvalues;
                                else
                                {
                                    BGC::BoundaryFluidElement mybgc;
                                    mybgc.rho = myvalues.alpha1rho1+myvalues.alpha2rho2;
                                    // IMPORTANT REMARK: if we need the dummy for something else
                                    // we can recompute the pressure from energy
                                    mybgc.pressure = myvalues.dummy;
                                    const double myrhoinv = 1.0/mybgc.rho;
                                    assert(!isnan(mybgc.rho));
                                    assert(!isnan(mybgc.pressure));
                                    assert(!isnan(myrhoinv));

                                    if (iz<0 && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY);
                                        delta_ix = 0;
                                        delta_iy = 0;
                                        delta_iz = 3;
                                        bnx = TBlock::sizeX;
                                        bny = TBlock::sizeY;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir2_side0[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir2_side0[bgb_index].fe[bgc_index].u=myvalues.rw*myrhoinv;
                                    }
                                    else if  (iz>=TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY);
                                        delta_ix = 0;
                                        delta_iy = 0;
                                        delta_iz = -TBlock::sizeZ;
                                        bnx = TBlock::sizeX;
                                        bny = TBlock::sizeY;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir2_side1[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir2_side1[bgb_index].fe[bgc_index].u=myvalues.rw*myrhoinv;
                                    }
                                    else if (iz>=0 && iz<TBlock::sizeZ && iy<0 && ix>=0 && ix<TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ);
                                        delta_iy = 3;
                                        delta_ix = 0;
                                        delta_iz = 0;
                                        bnx = TBlock::sizeX;
                                        bny = 3;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir1_side0[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir1_side0[bgb_index].fe[bgc_index].u=myvalues.rv*myrhoinv;
                                    }
                                    else if (iz>=0 && iz<TBlock::sizeZ && iy>=TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ);
                                        delta_iy = -TBlock::sizeY;
                                        delta_ix = 0;
                                        delta_iz = 0;
                                        bnx = TBlock::sizeX;
                                        bny = 3;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir1_side1[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir1_side1[bgb_index].fe[bgc_index].u=myvalues.rv*myrhoinv;
                                    }
                                    else if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix<0)
                                    {
                                        const int bgb_index = info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ);
                                        delta_ix = 3;
                                        delta_iy = 0;
                                        delta_iz = 0;
                                        bnx = 3;
                                        bny = TBlock::sizeY;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir0_side0[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir0_side0[bgb_index].fe[bgc_index].u=myvalues.ru*myrhoinv;
                                    }
                                    else if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ);
                                        delta_ix = -TBlock::sizeX;
                                        delta_iy = 0;
                                        delta_iz = 0;
                                        bny = TBlock::sizeY;
                                        bnx = 3;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir0_side1[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir0_side1[bgb_index].fe[bgc_index].u=myvalues.ru*myrhoinv;
                                    }
                                }
#endif
                            }
                }
            }
        }
    }

    if (m_clData.laplacep)
    {
        if (this->isroot)
            std::cout << "pressure relaxation for initial field currenrly not supported: ask FABIAN for further details" << std::endl;
        // _relax();
        // _set_energy(*grid);
    }

    if (this->isroot)
        cout << "done." << endl;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_Cloud<TGrid,TStepper,TSlice>::_parameter_cloud()
{
    // this data is defined in TestLabs.h
    CloudData::rho1 = this->parser("rho1").asDouble(1000.0);
    CloudData::rho2 = this->parser("rho2").asDouble(1.0);
    CloudData::p1   = this->parser("p1").asDouble(100.0);
    CloudData::p2   = this->parser("p2").asDouble(0.0234);

    // pack relevant cloud data passed to pressure and alpha handles
    m_clData.R0      = this->parser("-R0").asDouble(1.0);
    m_clData.extent[0]  = Simulation_Environment::extents[0];
    m_clData.extent[1]  = Simulation_Environment::extents[1];
    m_clData.extent[2]  = Simulation_Environment::extents[2];
    m_clData.epsilon = this->parser("-eps-sharp").asDouble(0.75);
    if (m_clData.epsilon <= 0.0)
        m_clData.epsilon = 0.75;
    m_clData.buff     = this->parser("-buffer").asDouble(0.0);
    m_clData.laplacep = this->parser("laplacep").asBool(false);
    m_clData.integral = this->parser("-ic-integral").asBool(true);

    // assign pressure and alpha handles
    const std::string ic_p = this->parser("-ic-pressure-cloud").asString("eqcloud");
    const bool cosine = this->parser("-cosine").asBool(true);
    pressure_handle = NULL;
    if (ic_p == "jump")
        pressure_handle = &_pressure_jump;
    else if (ic_p == "tiwari")
        pressure_handle = &_pressure_tiwari;
    else if (ic_p == "eqcloud")
    {
        if (cosine) pressure_handle = &_pressure_eqCloud<1>; // cos
        else        pressure_handle = &_pressure_eqCloud<0>; // tanh
        this->parser.set_strict_mode();
        m_clData.prox = this->parser("-proximity").asDouble();
        this->parser.unset_strict_mode();
    }
    else if (ic_p == "tiwari-cloud")
    {
        if (cosine) pressure_handle = &_pressure_tiwariCloud<1>; // cos
        else        pressure_handle = &_pressure_tiwariCloud<0>; // tanh
        this->parser.set_strict_mode();
        m_clData.prox = this->parser("-proximity").asDouble();
        this->parser.unset_strict_mode();
    }
    assert(pressure_handle != NULL);

    alpha_handle = NULL;
    if (ic_p != "tiwari")
        alpha_handle = &_alpha_std;
    else
        alpha_handle = &_alpha_tiwari;
    assert(alpha_handle != NULL);
}

#endif /* TEST_CLOUD_H_XEYR5OG6 */
