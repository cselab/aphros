/*
 *  Test_SICCloud.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger on 4/24/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_SICCLOUD_H_SH6ARFKA
#define TEST_SICCLOUD_H_SH6ARFKA

#include <omp.h>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

#include "Test_SteadyState.h"
#include "Shapes.h"

template <typename T>
class ShockType;


template <typename TGrid, typename TStepper, template <typename> class TSlice=Slice>
class Test_SICCloud: public virtual Test_SteadyState<TGrid,TStepper,TSlice>
{
public:
    Test_SICCloud(ArgumentParser& parser) :
        Test_SteadyState<TGrid,TStepper,TSlice>(parser), m_shock(NULL) { }
    virtual ~Test_SICCloud()
    {
        if (m_shock)
            delete m_shock;
    }


protected:
    typedef typename TGrid::BlockType TBlock;
    typedef typename TBlock::ElementType TElement;

    Real m_epsilon, m_mach, m_pressureRatio, m_shockSpeed;
    ShockType<TElement>* m_shock;

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////   TEST SHOCK INDUCED COLLAPSE   ///////////////\n");
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
        this->_parameter_sic_cloud();
    }

    virtual void _init()
    {
        Test_SteadyState<TGrid,TStepper,TSlice>::_init();
        this->_init_sic_cloud();
    }

    virtual void _ic();

    virtual void _analysis()
    {
        if (this->ANALYSISPERIOD != 0 && this->step_id%this->ANALYSISPERIOD == 0 && !this->bRESTART)
        {
            this->profiler.push_start("ANALYSIS");
            this->_sic_cloud_simple_analysis();
            this->profiler.pop_stop();
        }
    }

    void _parameter_sic_cloud();

    void _init_sic_cloud() const
    {
        if (this->isroot)
        {
            std::cout.setf(std::ios::scientific, std::ios::floatfield);
            std::cout << "SIC Cloud Initial condition..." << std::endl;
            std::cout << "INITIAL SHOCK" << std::endl;
            std::cout << '\t' << "p-Ratio         = " << m_pressureRatio << std::endl;
            std::cout << '\t' << "Mach            = " << m_mach << std::endl;
            std::cout << '\t' << "Shock speed     = " << m_shockSpeed << std::endl;
            std::cout << '\t' << "Shock direction = (" << m_shock->nx() << ", " << m_shock->ny() << ", " << m_shock->nz() << ")" << std::endl;
            std::cout << '\t' << "Point on shock  = (" << m_shock->Sx() << ", " << m_shock->Sy() << ", " << m_shock->Sz() << ")" << std::endl << std::endl;
            std::cout << "INITIAL MATERIALS" << std::endl;
            std::cout << '\t' << "MATERIAL 0:" << std::endl;
            std::cout << '\t' << "Gamma = " << Simulation_Environment::GAMMA1 << std::endl;
            std::cout << '\t' << "P_c   = " << Simulation_Environment::PC1 << std::endl;
            std::cout << "\t\t" << "Pre-Shock:" << std::endl;
            std::cout << "\t\t\t" << "rho = " << ShockType<TElement>::preRho << std::endl;
            std::cout << "\t\t\t" << "u   = " << ShockType<TElement>::preU << std::endl;
            std::cout << "\t\t\t" << "p   = " << ShockType<TElement>::preP << std::endl;
            std::cout << "\t\t" << "Post-Shock:" << std::endl;
            std::cout << "\t\t\t" << "rho = " << ShockType<TElement>::postRho << std::endl;
            std::cout << "\t\t\t" << "u   = " << ShockType<TElement>::postU << std::endl;
            std::cout << "\t\t\t" << "p   = " << ShockType<TElement>::postP << std::endl;
            std::cout << '\t' << "MATERIAL 1:" << std::endl;
            std::cout << '\t' << "Gamma = " << Simulation_Environment::GAMMA2 << std::endl;
            std::cout << '\t' << "P_c   = " << Simulation_Environment::PC2 << std::endl;
            std::cout << "\t\t" << "Pre-Shock:" << std::endl;
            std::cout << "\t\t\t" << "rho = " << ShockType<TElement>::bubbleRho << std::endl;
            std::cout << "\t\t\t" << "u   = " << ShockType<TElement>::bubbleU << std::endl;
            std::cout << "\t\t\t" << "p   = " << ShockType<TElement>::bubbleP << std::endl << std::endl;
        }
    }

    void _sic_cloud_simple_analysis();
    void _moving_shock_default();
    void _stationary_shock_rpp();
    void _stationary_shock_upp();
};


///////////////////////////////////////////////////////////////////////////////
// Shock Types
///////////////////////////////////////////////////////////////////////////////
template <typename TElement>
class ShockType
{
    virtual double _is_shock(const double pos[3]) const = 0;

protected:
    double m_nx, m_ny, m_nz; // shock normal
    double m_Sx, m_Sy, m_Sz; // point on shock surface

public:
    static double preRho, preU, preP;
    static double postRho, postU, postP;
    static double bubbleRho, bubbleU, bubbleP;

    ShockType(const double n[3], const double S[3]) :
        m_nx(n[0]), m_ny(n[1]), m_nz(n[2]), m_Sx(S[0]), m_Sy(S[1]), m_Sz(S[2]) { }

    virtual TElement set_shockBubble(const double pos[3], const double bubble = 0.0) const = 0;
    inline double nx() const { return m_nx; }
    inline double ny() const { return m_ny; }
    inline double nz() const { return m_nz; }
    inline double Sx() const { return m_Sx; }
    inline double Sy() const { return m_Sy; }
    inline double Sz() const { return m_Sz; }
};


template <typename T> double ShockType<T>::preRho = -1;
template <typename T> double ShockType<T>::preU   =  0;
template <typename T> double ShockType<T>::preP   = -1;

template <typename T> double ShockType<T>::postRho = -1;
template <typename T> double ShockType<T>::postU   =  0;
template <typename T> double ShockType<T>::postP   = -1;

template <typename T> double ShockType<T>::bubbleRho = -1;
template <typename T> double ShockType<T>::bubbleU   =  0;
template <typename T> double ShockType<T>::bubbleP   = -1;


// normal shock
template <typename TElement>
class NormalShock : public ShockType<TElement>
{
    virtual double _is_shock(const double pos[3]) const
    {
        // shock normal points into pre-shock region (naturally)
        const double d = this->m_nx*(pos[0] - this->m_Sx)
                + this->m_ny*(pos[1] - this->m_Sy)
                + this->m_nz*(pos[2] - this->m_Sz);

        return Simulation_Environment::heaviside(d);
    }

    const double m_epsilon;

public:
    NormalShock(const double n[3], const double S[3], const double eps=0) : ShockType<TElement>(n,S), m_epsilon(eps)
    {
        // normalize shock normal vector
        const double mag = std::sqrt(std::pow(this->m_nx, 2) + std::pow(this->m_ny, 2) + std::pow(this->m_nz, 2));
        assert(mag > 0);
        this->m_nx /= mag;
        this->m_ny /= mag;
        this->m_nz /= mag;
    }

    virtual TElement set_shockBubble(const double pos[3], const double bubble) const
    {
        const double G1 = 1.0/(Simulation_Environment::GAMMA1-1.0);
        const double G2 = 1.0/(Simulation_Environment::GAMMA2-1.0);
        const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
        const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;

        TElement out;
        const double shock = this->_is_shock(pos);

        // WARNING: THIS ONLY WORKS IF THE NORMAL VECTOR OF THE PLANAR SHOCK IS
        // colinear WITH EITHER UNIT VECTOR OF THE GLOBAL COORDINATE SYSTEM!
        const double ar1 = shock*this->postRho + (1-shock)*this->preRho*(1-bubble);
        const double ar2 = (1-shock)*this->bubbleRho*bubble;
        const double rhoMix = ar1 + ar2;
        out.alpha1rho1 = ar1;
        out.alpha2rho2 = ar2;

        const double u = shock*this->postU + (1-shock)*(this->bubbleU*bubble + this->preU*(1-bubble));
        out.ru = u*rhoMix;
        out.rv = 0.0;
        out.rw = 0.0;

        const double a2 = std::min(1.0-m_epsilon, std::max(m_epsilon, bubble));
        const double a1 = 1.0 - a2;
        const double pressure  = shock*this->postP + (1-shock)*(this->bubbleP*bubble + this->preP*(1-bubble));
        const double gmixm1Inv = a1*G1 + a2*G2;
        const double pcmix     = a1*G1*F1 + a2*G2*F2;
        const double Ekin      = 0.5*rhoMix*u*u;
        out.energy = gmixm1Inv*pressure + pcmix + Ekin;
        out.alpha2 = a2;
        out.dummy = 0.0;

        return out;
    }
};

// spherical shock
template <typename TElement>
class SphericalShock : public ShockType<TElement>
{
    virtual double _is_shock(const double pos[3]) const
    {
        const double rx = (pos[0] - this->m_Sx);
        const double ry = (pos[1] - this->m_Sy);
        const double rz = (pos[2] - this->m_Sz);
        const double r  = std::sqrt(rx*rx + ry*ry + rz*rz);

        // TODO: (fabianw@mavt.ethz.ch; Fri 01 Apr 2016 12:21:25 PM CEST) check
        // how smooth shock interface behaves
        return Simulation_Environment::heaviside(m_radius - r);
        /* return Simulation_Environment::heaviside_smooth(m_radius - r); */
    }

    const double m_radius;
    const double m_epsilon;

public:
    SphericalShock(const double n[3], const double S[3], const double r, const double eps=0) : ShockType<TElement>(n,S), m_radius(r), m_epsilon(eps)
    { }

    virtual TElement set_shockBubble(const double pos[3], const double bubble) const
    {
        const double G1 = 1.0/(Simulation_Environment::GAMMA1-1.0);
        const double G2 = 1.0/(Simulation_Environment::GAMMA2-1.0);
        const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
        const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;

        TElement out;
        const double shock = this->_is_shock(pos);

        // density
        const double ar1 = shock*this->postRho + (1-shock)*this->preRho*(1-bubble);
        const double ar2 = (1-shock)*this->bubbleRho*bubble;
        const double rhoMix = ar1 + ar2;
        out.alpha1rho1 = ar1;
        out.alpha2rho2 = ar2;

        // momentum
        const double rx = (pos[0] - this->m_Sx);
        const double ry = (pos[1] - this->m_Sy);
        const double rz = (pos[2] - this->m_Sz);
        const double r  = std::sqrt(rx*rx + ry*ry + rz*rz);
        double theta, phi, ux, uy, uz;
        if (r > 0.0)
        {
            // we treat postU as the velocity magnitude observed in the post
            // shock region
            theta = std::atan2(ry, rx);
            phi   = std::acos(rz/r);
            ux    = -this->postU*sin(phi)*cos(theta);
            uy    = -this->postU*sin(phi)*sin(theta);
            uz    = -this->postU*cos(phi);
        }
        else
        {
            // what is happening at the singular point?  in our case, we assume
            // a quiescent fluid, so we set ux = uy = uz = 0
            ux = uy = uz = 0.0;
        }

        // TODO: (fabianw@mavt.ethz.ch; Fri 01 Apr 2016 12:27:57 PM CEST) this
        // assumes zero momentum in the pre-shock region!
        const double u = (shock*ux + (1-shock)*(this->bubbleU*bubble + this->preU*(1-bubble)));
        const double v = (shock*uy + (1-shock)*(this->bubbleU*bubble + this->preU*(1-bubble)));
        const double w = (shock*uz + (1-shock)*(this->bubbleU*bubble + this->preU*(1-bubble)));
        out.ru = u*rhoMix;
        out.rv = v*rhoMix;
        out.rw = w*rhoMix;

        // Energia
        const double a2 = std::min(1.0-m_epsilon, std::max(m_epsilon, bubble));
        const double a1 = 1.0 - a2;
        const double pressure  = shock*this->postP + (1-shock)*(this->bubbleP*bubble + this->preP*(1-bubble));
        const double gmixm1Inv = a1*G1 + a2*G2;
        const double pcmix     = a1*G1*F1 + a2*G2*F2;
        const double Ekin      = 0.5*rhoMix*(u*u + v+v + w*w);
        out.energy = gmixm1Inv*pressure + pcmix + Ekin;
        out.alpha2 = a2;
        out.dummy = 0.0;

        return out;
    }
};


///////////////////////////////////////////////////////////////////////////////
// Helper
///////////////////////////////////////////////////////////////////////////////
template <typename T>
static inline T get_IC(const double pos[3], const ShockType<T> * const shock, const std::vector<shape> * const blockShapes = NULL)
{
    T IC;

    /* *
     * For the IC, bubbles are not allowed to be located in the post shock
     * region.  No check is being performed.
     * */
    if (!blockShapes)
    {
        // If blockShapes == NULL, it is assumed there are no shapes in this
        // block.  Therefore, only the shock wave is treated for initial
        // conditions.
        IC = shock->set_shockBubble(pos);
    }
    else
    {
        // There are shapes within this block to be taken care of.
        const double bubble = eval(*blockShapes, pos);
        IC = shock->set_shockBubble(pos, bubble);
    }

    return IC;
}


template<typename T>
static inline T integral(const double p[3], const double h, const ShockType<T> * const myshock, const std::vector<shape> * const blockShapes = NULL) // h should be half the grid spacing h
{
    T samples[3][3][3];
    T zintegrals[3][3];
    T yzintegrals[3];

    const double x0[3] = {p[0] - h, p[1] - h, p[2] - h};

    for(int iz=0; iz<3; iz++)
        for(int iy=0; iy<3; iy++)
            for(int ix=0; ix<3; ix++)
            {
                const double mypos[3] = {x0[0]+ix*h, x0[1]+iy*h, x0[2]+iz*h};
                samples[iz][iy][ix] = get_IC<T>(mypos, myshock, blockShapes);
            }

    for(int iy=0; iy<3; iy++)
        for(int ix=0; ix<3; ix++)
            zintegrals[iy][ix] = (1/6.) * samples[0][iy][ix]+(2./3) * samples[1][iy][ix]+(1/6.)* samples[2][iy][ix];

    for(int ix=0; ix<3; ix++)
        yzintegrals[ix] = (1./6)*zintegrals[0][ix] + (2./3)*zintegrals[1][ix]+(1./6)*zintegrals[2][ix];

    return (1./6) * yzintegrals[0]+(2./3) * yzintegrals[1]+(1./6)* yzintegrals[2];
}


///////////////////////////////////////////////////////////////////////////////
// Class Implementation
///////////////////////////////////////////////////////////////////////////////
template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SICCloud<TGrid,TStepper,TSlice>::_ic()
{
    // set up bubbles
    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();
#ifdef _NONUNIFORM_BLOCK_
    const double h = this->m_nonuniform->minimum_cell_width();
#else
    const double h = vInfo[0].h_gridpoint;
#endif /* _NONUNIFORM_BLOCK_ */
    std::string b_file = this->parser("-cloud").asString("cloud.dat");
    Seed<shape> myseed;
    myseed.make_shapes(b_file, h, this->isroot);

    // set fields
#pragma omp parallel
    {
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif

        if (myseed.get_shapes().size() == 0)
        {
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
                            b(ix,iy,iz) = get_IC<TElement>(pos, m_shock);
                            // b(ix,iy,iz) = integral<TElement>(pos, 0.5*h);
                        }
            }
        }
        else
#pragma omp for
            for(int i=0; i<(int)vInfo.size(); i++)
            {
                BlockInfo info = vInfo[i];
                double myextent[3] = { info.block_extent[0], info.block_extent[1], info.block_extent[2] } ;
                TBlock& b = *(TBlock*)info.ptrBlock;
                std::vector<shape> myshapes = myseed.get_shapes().size() ? myseed.retain_shapes(info.origin, myextent).get_shapes() : myseed.get_shapes();
                const bool isempty = myshapes.size() == 0;

                if (isempty)
                {
                    for(int iz=0; iz<TBlock::sizeZ; iz++)
                        for(int iy=0; iy<TBlock::sizeY; iy++)
                            for(int ix=0; ix<TBlock::sizeX; ix++)
                            {
                                double pos[3];
                                info.pos(pos, ix, iy, iz);
                                b(ix,iy,iz) = get_IC<TElement>(pos, m_shock);
                                // b(ix,iy,iz) = integral<TElement>(pos, 0.5*h);
                            }
                }
                else
                {
                    for(int iz=0; iz<TBlock::sizeZ; iz++)
                        for(int iy=0; iy<TBlock::sizeY; iy++)
                            for(int ix=0; ix<TBlock::sizeX; ix++)
                            {
                                double pos[3];
                                info.pos(pos, ix, iy, iz);
                                // b(ix,iy,iz) = get_IC<TElement>(pos, &myshapes);
                                b(ix,iy,iz) = integral<TElement>(pos, 0.5*h, m_shock, &myshapes);
                            }
                }
            }
    }
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SICCloud<TGrid,TStepper,TSlice>::_parameter_sic_cloud()
{
    // mixture epsilon
    m_epsilon = this->parser("-epsilon").asDouble(0.0);

    // bubbles (in pre-shock region)
    ShockType<TElement>::bubbleRho = this->parser("-rhoB").asDouble(1);
    ShockType<TElement>::bubbleU   = this->parser("-uB").asDouble(0);
    ShockType<TElement>::bubbleP   = this->parser("-pB").asDouble(1);

    // setup shock
    const std::string shocktype(this->parser("shocktype").asString("normal"));

    const bool bSrpp = this->parser("stationaryshock_rpp").asBool(false);
    const bool bSupp = this->parser("stationaryshock_upp").asBool(false);

    double n[3] = {0};
    double S[3] = {0};
    if (shocktype == "normal")
    {
        // shock orienation
        n[0] = this->parser("-shockNx").asDouble(1.0);
        n[1] = this->parser("-shockNy").asDouble(0.0);
        n[2] = this->parser("-shockNz").asDouble(0.0);
        // point on shock
        S[0] = this->parser("-shockSx").asDouble(0.05);
        S[1] = this->parser("-shockSy").asDouble(0.0);
        S[2] = this->parser("-shockSz").asDouble(0.0);
        m_shock = new NormalShock<TElement>(n,S,m_epsilon);
    }
    else if (shocktype == "spherical")
    {
        // initial shock radius
        const double Rs = this->parser("-shockradius").asDouble(0.4);
        // singular point (assume domain center)
        S[0] = this->parser("-shockSx").asDouble(0.5*Simulation_Environment::extents[0]);
        S[1] = this->parser("-shockSy").asDouble(0.5*Simulation_Environment::extents[1]);
        S[2] = this->parser("-shockSz").asDouble(0.5*Simulation_Environment::extents[2]);
        m_shock = new SphericalShock<TElement>(n,S,Rs,m_epsilon);
    }
    else
    {
        std::cerr << "ERROR: unknown shock type \"" << shocktype << "\"" << std::endl;
        abort();
    }

    // init shock
    if ( (bSrpp || bSupp) && shocktype == "normal")
    {
        if (bSrpp) // specifying both does not make sense, so rpp goes first
            _stationary_shock_rpp();
        else
            _stationary_shock_upp();
    }
    else
        this->_moving_shock_default();
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SICCloud<TGrid,TStepper,TSlice>::_sic_cloud_simple_analysis()
{
    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();
    double rInt=0., uInt=0., vInt=0., wInt=0., eInt=0., vol=0., ke=0., mach_max=-HUGE_VAL, p_max=-HUGE_VAL, p_min = HUGE_VAL;
    const double h = vInfo[0].h_gridpoint;
    const double h3 = h*h*h;

    for(int i=0; i<(int)vInfo.size(); i++)
    {
        BlockInfo info = vInfo[i];
        TBlock& b = *(TBlock*)info.ptrBlock;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    const double rhoMix = b(ix,iy,iz).alpha1rho1 + b(ix,iy,iz).alpha2rho2;
                    const double a2 = b(ix,iy,iz).alpha2;
                    const double a1 = 1.0 - a2;
                    rInt += rhoMix;
                    uInt += b(ix,iy,iz).ru;
                    vInt += b(ix,iy,iz).rv;
                    wInt += b(ix,iy,iz).rw;
                    eInt += b(ix,iy,iz).energy;
                    vol  += a2;
                    ke   += 0.5/rhoMix * (b(ix,iy,iz).ru*b(ix,iy,iz).ru+b(ix,iy,iz).rv*b(ix,iy,iz).rv+b(ix,iy,iz).rw*b(ix,iy,iz).rw);

                    const double p = Simulation_Environment::pressure_5eq(rhoMix, b(ix,iy,iz).energy, b(ix,iy,iz).ru, b(ix,iy,iz).rv, b(ix,iy,iz).rw, a2);
                    const double _g1 = Simulation_Environment::GAMMA1;
                    const double _g2 = Simulation_Environment::GAMMA2;
                    const double _pc1 = Simulation_Environment::PC1;
                    const double _pc2 = Simulation_Environment::PC2;
                    const double c = std::sqrt(_g1*_g2*(p+_pc1)*(p+_pc2) / ((a1*_g2*(p+_pc2) + a2*_g1*(p+_pc1))*rhoMix));
                    const double velmag = std::sqrt(b(ix,iy,iz).ru*b(ix,iy,iz).ru+b(ix,iy,iz).rv*b(ix,iy,iz).rv+b(ix,iy,iz).rw*b(ix,iy,iz).rw)/rhoMix;

                    mach_max = std::max(mach_max, velmag/c);
                    p_max = std::max(p_max, p);
                    p_min = std::min(p_min, p);
                }
    }

    FILE * f = fopen("integrals.dat", "a");
    fprintf(f, "%d %e %e %e %e %e %e %e %e %e %e %e %e %e\n", this->step_id, this->t, this->dt, rInt*h3, uInt*h3,
            vInt*h3, wInt*h3, eInt*h3, vol*h3, ke*h3, mach_max, p_max, p_min, std::pow(0.75*vol*h3/M_PI,1./3.));
    fclose(f);
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SICCloud<TGrid,TStepper,TSlice>::_moving_shock_default()
{
    const double gamma = Simulation_Environment::GAMMA1;
    const double pc    = Simulation_Environment::PC1;
    const double G1 = 1.0 / (gamma - 1.0);

    // liquid (in pre-shock region)
    ShockType<TElement>::preRho  = this->parser("-rho1").asDouble(1000);
    ShockType<TElement>::preU    = this->parser("-u1").asDouble(0);
    ShockType<TElement>::preP    = this->parser("-p1").asDouble(1);

    const bool bWithMach = this->parser.check("-mach");
    // pressure ratio over shock
    if (bWithMach)
        m_mach = this->parser("-mach").asDouble(1.027976);
    else
        m_pressureRatio = this->parser("-pressureratio").asDouble(400); // default

    /* *
     * Compute post shock states based on Appendix A of E. Johnsen "Numerical
     * Simulations of Non-Spherical Bubble Collapse", PhD Thesis, Caltech 2007.
     * */
    const double rho1  = ShockType<TElement>::preRho;
    const double u1    = ShockType<TElement>::preU;
    const double p1    = ShockType<TElement>::preP;
    if (bWithMach)
        m_pressureRatio = 1.0 + 2.0*gamma/(gamma+1.0)*(m_mach*m_mach - 1.0)*(1.0 + pc/p1);
    const double p2    = p1*m_pressureRatio;
    const double tmp1  = (gamma + 1.0)/(gamma - 1.0);
    const double tmp2  = (p2 + pc)/(p1 + pc);

    const double c1 = std::sqrt(gamma*(p1 + pc)/rho1); // pure liquid
    if (!bWithMach)
        m_mach = std::sqrt((gamma+1.0)/(2.0*gamma)*(m_pressureRatio - 1.0)*p1/(p1 + pc) + 1.0);

    const double u2   = u1 + c1*(m_pressureRatio - 1.0)*p1/(gamma*(p1 + pc)*m_mach);
    const double rho2 = rho1*(tmp1*tmp2 + 1.0)/(tmp1 + tmp2);
    m_shockSpeed = u1 + c1*m_mach;

    const std::string shocktype(this->parser("shocktype").asString("normal"));
    double sign = 1.0;
    if (shocktype == "normal")
    {
        const double nx = this->m_shock->nx();
        assert(nx != 0.0);
        sign = ((0.0 < nx) - (nx < 0.0)) * sign;
    }

    ShockType<TElement>::postRho = rho2;
    ShockType<TElement>::postU   = sign*u2;
    ShockType<TElement>::postP   = p2;

    // assign values at boundaries (post-shock conditions)
    SICCloudBCData::boundaryElement0.alpha1rho1 = rho2;
    SICCloudBCData::boundaryElement0.alpha2rho2 = 0.0;
    SICCloudBCData::boundaryElement0.ru         = sign*rho2*u2;
    SICCloudBCData::boundaryElement0.rv         = 0.0;
    SICCloudBCData::boundaryElement0.rw         = 0.0;
    SICCloudBCData::boundaryElement0.energy = p2*G1 + pc*gamma*G1 + 0.5*rho2*u2*u2; // v = w = 0 !
    SICCloudBCData::boundaryElement0.alpha2 = m_epsilon;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SICCloud<TGrid,TStepper,TSlice>::_stationary_shock_rpp()
{
    const double gamma = Simulation_Environment::GAMMA1;
    const double pc    = Simulation_Environment::PC1;
    const double G1 = 1.0 / (gamma - 1.0);

    // specify pre-shock rho and pressure, as well as post-shock pressure
    ShockType<TElement>::preRho  = this->parser("-rho1").asDouble(1);
    ShockType<TElement>::preP    = this->parser("-p1").asDouble(1);
    ShockType<TElement>::postP   = this->parser("-p2").asDouble(100);

    const double rho1 = ShockType<TElement>::preRho;
    const double p1   = ShockType<TElement>::preP;
    const double p2   = ShockType<TElement>::postP;

    const double u1   = std::sqrt( 1.0/(rho1*G1)*(p2*G1 + gamma*pc*G1 + 0.5*(p1 + p2)) );
    const double u2   = 1.0/(rho1*u1)*(rho1*u1*u1 + p1 - p2);
    const double rho2 = rho1*u1/u2;

    m_pressureRatio = p2/p1;
    const double c1 = std::sqrt(gamma*(p1 + pc)/rho1); // pure liquid here
    m_mach = std::sqrt((gamma+1.0)/(2.0*gamma)*(m_pressureRatio - 1.0)*p1/(p1 + pc) + 1.0);
    m_shockSpeed = u1 + c1*m_mach;

    ShockType<TElement>::bubbleU = u1;
    ShockType<TElement>::preU    = u1;
    ShockType<TElement>::postRho = rho2;
    ShockType<TElement>::postU   = u2;

    // assign values at boundaries (pre-shock conditions)
    SICCloudBCData::boundaryElement0.alpha1rho1 = rho1;
    SICCloudBCData::boundaryElement0.alpha2rho2 = 0.0;
    SICCloudBCData::boundaryElement0.ru         = rho1*u1;
    SICCloudBCData::boundaryElement0.rv         = 0.0;
    SICCloudBCData::boundaryElement0.rw         = 0.0;
    SICCloudBCData::boundaryElement0.energy = p1*G1 + pc*gamma*G1 + 0.5*rho1*u1*u1; // v = w = 0 !
    SICCloudBCData::boundaryElement0.alpha2 = m_epsilon;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SICCloud<TGrid,TStepper,TSlice>::_stationary_shock_upp()
{
    const double gamma = Simulation_Environment::GAMMA1;
    const double pc    = Simulation_Environment::PC1;
    const double G1 = 1.0 / (gamma - 1.0);

    const double densityRatio = this->parser("-densityratio").asDouble(1000.0); // rho1/rhoB

    // specify pre-shock u and pressure, as well as post-shock pressure
    ShockType<TElement>::preU  = this->parser("-u1").asDouble(1);
    ShockType<TElement>::preP  = this->parser("-p1").asDouble(1);
    ShockType<TElement>::postP = this->parser("-p2").asDouble(100);

    const double u1 = ShockType<TElement>::preU;
    const double p1 = ShockType<TElement>::preP;
    const double p2 = ShockType<TElement>::postP;

    const double cf = gamma*pc*G1 + 0.5*(p1 + p2);
    const double u2 = u1*(p1*G1 + cf) / (p2*G1 + cf);
    const double rho2 = (p2 - p1) / (u2*(u1 - u2));
    const double rho1 = rho2*u2/u1;

    m_pressureRatio = p2/p1;
    const double c1 = std::sqrt(gamma*(p1 + pc)/rho1); // pure liquid here
    m_mach = std::sqrt((gamma+1.0)/(2.0*gamma)*(m_pressureRatio - 1.0)*p1/(p1 + pc) + 1.0);
    m_shockSpeed = u1 + c1*m_mach;

    ShockType<TElement>::bubbleU = u1;
    ShockType<TElement>::bubbleRho = rho1 / densityRatio;

    ShockType<TElement>::preRho  = rho1;
    ShockType<TElement>::postRho = rho2;
    ShockType<TElement>::postU   = u2;

    // assign values at boundaries (pre-shock conditions)
    SICCloudBCData::boundaryElement0.alpha1rho1 = rho1;
    SICCloudBCData::boundaryElement0.alpha2rho2 = 0.0;
    SICCloudBCData::boundaryElement0.ru         = rho1*u1;
    SICCloudBCData::boundaryElement0.rv         = 0.0;
    SICCloudBCData::boundaryElement0.rw         = 0.0;
    SICCloudBCData::boundaryElement0.energy = p1*G1 + pc*gamma*G1 + 0.5*rho1*u1*u1; // v = w = 0 !
    SICCloudBCData::boundaryElement0.alpha2 = m_epsilon;
}
#endif /* TEST_SICCLOUD_H_SH6ARFKA */
