/*
 *  FlowStep_LSRK3.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <StencilInfo.h>

#include <Source_CPP.h>
#include <Update.h>
#include <Update_State_CPP_5eq.h>
#include <MaxSpeedOfSound_CPP_5eq.h>
#include <Convection_CPP_HLLC_5eq.h>
#include <Convection_CPP_HLLC_5eq_nonuniform.h>
#include <Diffusion_CPP.h>
#include <SurfaceTension_CPP.h>
#include <InterfaceSharpening_CPP_cell_compact_5eq.h>
#include <InterfaceSharpening_CPP_5eq.h>
#include <MaxInterfaceVel_CPP_5eq.h>

#if defined(_QPX_) || defined(_QPXEMU_)
#include <Update_QPX.h>
#include <MaxSpeedOfSound_QPX_5eq.h>
#include <Convection_QPX_HLLC_5eq.h>
#if defined(_BGQ_)
#else
//#include <Diffusion_QPX.h> // dummy for now
#endif
#endif
#ifdef _QPX_
//here are for cycle counting
#include <ucontext.h>
#include <signal.h>
#include <sys/time.h>
#include <errno.h>
#include "spi/include/upci/upci.h"
#endif

#ifdef _USE_HPM_
#include <mpi.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#else
#define HPM_Start(x)
#define HPM_Stop(x)
#endif

#include "BlockProcessor_OMP.h"
#include "Types.h"
#include "NonUniform.h"

/* TODO: (fabianw; Wed 03 Jun 2015 07:41:07 PM CEST) see inside Tests.h */
#include "Tests.h"


template < typename TIVEL>
Real _computeMaxIVel_OMP(Grid_t& grid)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    const int N = vInfo.size();
    const BlockInfo * const ary = &vInfo.front();

#ifdef _USE_HPM_
    HPM_Start("dt");
#endif

#if (_OPENMP < 201107)
    Real * tmp = NULL;
    int error = posix_memalign((void**)&tmp, std::max(8, _ALIGNBYTES_), sizeof(Real) * N);
    assert(error == 0);

    Real * const local_ivel = tmp;

#pragma omp parallel
    {

#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif

        TIVEL kernel;
#pragma omp for schedule(runtime)
        for (int i=0; i<N; ++i)
        {
            Block_t & block = *(Block_t *)ary[i].ptrBlock;
            local_ivel[i] =  kernel.compute(&block.data[0][0][0].alpha1rho1, Block_t::gptfloats); //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
        }
    }

    Real global_ivel = local_ivel[0];

#pragma omp parallel
    {
        Real mymax = local_ivel[0];

#pragma omp for schedule(runtime)
        for(int i=0; i<N; ++i)
            mymax = max(local_ivel[i], mymax);

#pragma omp critical
        {
            global_ivel = max(global_ivel, mymax);
        }
    }

    free(tmp);

#else   /* FUSED_MAXSOS */

    Real global_ivel = 0;

#pragma omp parallel
    {
        TIVEL kernel;

#pragma omp for schedule(runtime) reduction(max:global_ivel)
        for (size_t i=0; i<N; ++i)
        {
            Block_t & block = *(Block_t *)ary[i].ptrBlock;
            global_ivel =  max(global_ivel, kernel.compute(&block.data[0][0][0].alpha1rho1, Block_t::gptfloats)); //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
        }
    }

#endif

#ifdef _USE_HPM_
    HPM_Stop("dt");
#endif

    return global_ivel;
}


namespace LSRK3data
{
    extern Real smoothlength ;
    extern int verbosity;
    extern float PEAKPERF_CORE, PEAKBAND;
    extern int NCORES;
    extern string dispatcher;
    extern int ReportFreq;
    extern int step_id;

    extern int sponge, sponge_width[3];
    extern Real sponge_pref;

    extern int state;
    extern Real state_alpha, state_beta;

    template < typename Kernel , typename Lab>
    struct FlowStep
    {
        StencilInfo stencil;

        Real a, dtinvh;

        int stencil_start[3];
        int stencil_end[3];

        FlowStep(Real a, Real dtinvh): a(a), dtinvh(dtinvh), stencil(-3,-3,-3,4,4,4, false, 7, 0,1,2,3,4,5,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
        }

        FlowStep(const FlowStep& c): a(c.a), dtinvh(c.dtinvh), stencil(-3,-3,-3,4,4,4, false, 7, 0,1,2,3,4,5,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
        }

        inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
        {
            //TODO: URSULA temporary extenxtion of constructor arguments
            Kernel kernel(a, dtinvh, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, Simulation_Environment::PC1, Simulation_Environment::PC2);

            const Real * const srcfirst = &lab(-3,-3,-3).alpha1rho1; //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
            const int labSizeRow = lab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
            Real * const destfirst = &o.tmp[0][0][0][0];

            kernel.compute(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
                           destfirst, Block_t::gptfloats, Block_t::sizeX, Block_t::sizeX*Block_t::sizeY);
#ifdef _QUANTIFY_DELTAP_
            const size_t blockID = info.blockID;
            const Real maxVal = kernel.getMaxPCorrection();
            const size_t cnt  = kernel.getNumberOfPCorrections();
            kernel.resetPCorrection();
            if (cnt)
            {
#pragma omp critical
                printf("QUANTIFY DELTAP: block=%d, count=%d, maxVal=%e\n", blockID, cnt, maxVal);
            }
#endif /* _QUANTIFY_DELTAP_ */
        }
    };

#ifdef _NONUNIFORM_BLOCK_
    template < typename Kernel , typename Lab>
    struct FlowStepNonuniform
    {
        StencilInfo stencil;

        Real a, dt;

        int stencil_start[3];
        int stencil_end[3];

        FlowStepNonuniform(Real a, Real dt): a(a), dt(dt), stencil(-3,-3,-3,4,4,4, false, 7, 0,1,2,3,4,5,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
        }

        FlowStepNonuniform(const FlowStepNonuniform& c): a(c.a), dt(c.dt), stencil(-3,-3,-3,4,4,4, false, 7, 0,1,2,3,4,5,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
        }

        inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
        {
            Kernel kernel(a, dt, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, Simulation_Environment::PC1, Simulation_Environment::PC2);

            const Real * const srcfirst = &lab(-3,-3,-3).alpha1rho1;
            const int labSizeRow = lab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
            Real * const destfirst = &o.tmp[0][0][0][0];

            kernel.compute_nonuniform(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
                    destfirst, Block_t::gptfloats, Block_t::sizeX, Block_t::sizeX*Block_t::sizeY,
                    o.coeffs_x[1], o.coeffs_x[0], o.coeffs_y[1], o.coeffs_y[0], o.coeffs_z[1], o.coeffs_z[0],
                    &o.invh_x[0], &o.invh_y[0], &o.invh_z[0]);
        }
    };
#endif /*  */


    template < typename Kernel>
    struct Source
    {
        Real a, dt;
        Real g[3], f[3];
        BlockInfo * ary;

        Source (Real a, Real dt, const Real gaccel[3], const Real volforce[3], BlockInfo * ary) :
            a(a), dt(dt), ary(ary)
        {
            g[0]=gaccel[0];
            g[1]=gaccel[1];
            g[2]=gaccel[2];
            f[0]=volforce[0];
            f[1]=volforce[1];
            f[2]=volforce[2];
        }

        Source (const Source& c): a(c.a), dt(c.dt), g(c.g), f(c.f), ary(c.ary) {}

        void omp(const int N)
        {
#pragma omp parallel
            {
#ifdef _USE_NUMA_
                const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
                const int mynode = omp_get_thread_num() / cores_per_node;
                numa_run_on_node(mynode);
#endif
                Kernel kernel(a, dt, g, f);

#pragma omp for schedule(runtime)
                for(int r=0; r<N; ++r)
                {
                    Block_t & block = *(Block_t *)ary[r].ptrBlock;
                    kernel.compute(&block.data[0][0][0].alpha1rho1, block.gptfloats, &block.tmp[0][0][0][0], block.gptfloats);
                }
            }
        }

     };


    template <typename Kernel>
    struct OneWayAcousticSource
    {
        typedef typename Kernel::SourceParameter param_t;
        Real t, dt;
        param_t p;
        BlockInfo * ary;

        OneWayAcousticSource(const Real t, const Real dt, const param_t& p, BlockInfo * ary) : t(t), dt(dt), p(p), ary(ary) {}
        OneWayAcousticSource(const OneWayAcousticSource& c): t(c.t), dt(c.dt), p(c.p), ary(c.ary) {}

        void omp(const int N)
        {
#pragma omp parallel
            {
#ifdef _USE_NUMA_
                const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
                const int mynode = omp_get_thread_num() / cores_per_node;
                numa_run_on_node(mynode);
#endif
                Kernel kernel(t, dt, p);

#pragma omp for schedule(runtime)
                for(int r=0; r<N; ++r)
                {
                    Block_t & block = *(Block_t *)ary[r].ptrBlock;
                    kernel.compute(&block.data[0][0][0].alpha1rho1, block.gptfloats, &block.tmp[0][0][0][0], block.gptfloats, ary[r]);
                }
            }
        }

     };

    template < typename Kernel, typename Lab >
    struct Diffusion
    {
        StencilInfo stencil;
        Real mu1, mu2, dtinvh;
        int stencil_start[3];
        int stencil_end[3];

// #ifndef _LIQUID_
//         Diffusion(const Real _dtinvh, const Real _mu1=0, const Real _mu2=0): dtinvh(_dtinvh), mu1(_mu1), mu2(_mu2), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,5)
// #else
//         Diffusion(const Real _dtinvh, const Real _mu1=0, const Real _mu2=0): dtinvh(_dtinvh), mu1(_mu1), mu2(_mu2), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,7)
// #endif
        Diffusion(const Real _dtinvh, const Real _mu1=0, const Real _mu2=0): dtinvh(_dtinvh), mu1(_mu1), mu2(_mu2), stencil(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }

// #ifndef _LIQUID_
//         Diffusion(const Diffusion& c): dtinvh(c.dtinvh), mu1(c.mu1), mu2(c.mu2), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,5)
// #else
//         Diffusion(const Diffusion& c): dtinvh(c.dtinvh), mu1(c.mu1), mu2(c.mu2), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,7)
// #endif
        Diffusion(const Diffusion& c): dtinvh(c.dtinvh), mu1(c.mu1), mu2(c.mu2), stencil(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }

        inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
        {
            Kernel kernel(1, mu1, mu2, info.h_gridpoint, LSRK3data::smoothlength, dtinvh);

            const Real * const srcfirst = &lab(-1,-1,-1).alpha1rho1;
            const int labSizeRow = lab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
            Real * const destfirst =  &o.tmp[0][0][0][0];
            kernel.compute(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
                           destfirst, Block_t::gptfloats, Block_t::sizeX, Block_t::sizeX*Block_t::sizeY);
        }
    };


    template < typename Kernel, typename Lab >
    struct SurfaceTension
    {
        StencilInfo stencil;
        Real sigma, dtinvh;
        int stencil_start[3];
        int stencil_end[3];

//#ifndef _LIQUID_
//        SurfaceTension(const Real _dtinvh, const Real _sigma=0.0): dtinvh(_dtinvh), sigma(_sigma), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,5)
//#else
//        SurfaceTension(const Real _dtinvh, const Real _sigma=0.0): dtinvh(_dtinvh), sigma(_sigma), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,7)
//#endif
        SurfaceTension(const Real _dtinvh, const Real _sigma=0.0): dtinvh(_dtinvh), sigma(_sigma), stencil(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }

//#ifndef _LIQUID_
//        SurfaceTension(const SurfaceTension& c): dtinvh(c.dtinvh), sigma(c.sigma), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,5)
//#else
//        SurfaceTension(const SurfaceTension& c): dtinvh(c.dtinvh), sigma(c.sigma), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,7)
//#endif
        SurfaceTension(const SurfaceTension& c): dtinvh(c.dtinvh), sigma(c.sigma), stencil(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }

        inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
        {
            Kernel kernel(1.0, dtinvh, info.h_gridpoint, sigma);

            const Real * const srcfirst = &lab(-1,-1,-1).alpha1rho1;
            const int labSizeRow = lab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
            Real * const destfirst = &o.tmp[0][0][0][0];
            kernel.compute(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
                           destfirst, Block_t::gptfloats, Block_t::sizeX, Block_t::sizeX*Block_t::sizeY);
        }
    };

    template < typename Kernel, typename Lab >
    struct InterfaceSharpening
    {
      StencilInfo stencil;
      Real dtinvh;
      Real epsilon, U0;
      int stencil_start[3];
      int stencil_end[3];

      InterfaceSharpening(const Real _dtinvh, const Real _epsilon=0, const Real _U0=0): dtinvh(_dtinvh), epsilon(_epsilon), U0(_U0), stencil(-1,-1,-1,2,2,2, true, 7, 0,1,2,3,4,5,6)
      {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
      }

      inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
      {
/*
         Kernel kernel(1, dtinvh, info.h_gridpoint, epsilon, U0, (Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1),
                       LSRK3data::gamma1*LSRK3data::pc1/(LSRK3data::gamma1-1), LSRK3data::gamma2*LSRK3data::pc2/(LSRK3data::gamma2-1));
*/
         Kernel kernel(1, dtinvh, epsilon, U0, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2,
                       Simulation_Environment::PC1, Simulation_Environment::PC2);

         const Real * const srcfirst = &lab(-1,-1,-1).alpha1rho1; //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
         const int labSizeRow = lab.template getActualSize<0>();
         const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
         Real * const destfirst =  &o.tmp[0][0][0][0];
         kernel.compute(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
                        destfirst, Block_t::gptfloats, Block_t::sizeX, Block_t::sizeX*Block_t::sizeY);
      }
    };


    template < typename Lab >
    struct UpdateBC
    {
        StencilInfo stencil;
        Real dt, a, b;
        int stencil_start[3];
        int stencil_end[3];

        UpdateBC(const Real _dt, const Real _a, const Real _b): dt(_dt), a(_a), b(_b), stencil(-3,-3,-3,4,4,4, false, 7, 0,1,2,3,4,5,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
        }

        UpdateBC(const UpdateBC& c): dt(c.dt), a(c.a), b(c.b), stencil(-3,-3,-3,4,4,4, false, 7, 0,1,2,3,4,5,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
        }

        inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
        {
          lab.apply_bc_update(info, dt, a, b);
        }
    };


    template < typename Kernel >
    struct Update
    {
        Real b;
        BlockInfo * ary;

    public:

        Update(float b, BlockInfo * ary): b(b), ary(ary) { }
        Update(const Update& c): b(c.b), ary(c.ary) { }

        void omp(const int N)
        {
#pragma omp parallel
            {
#ifdef _USE_NUMA_
                const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
                const int mynode = omp_get_thread_num() / cores_per_node;
                numa_run_on_node(mynode);
#endif
                Kernel kernel(b);

#pragma omp for schedule(runtime)
                for(int r=0; r<N; ++r)
                {
                    Block_t & block = *(Block_t *)ary[r].ptrBlock;
                    kernel.compute(&block.tmp[0][0][0][0], &block.data[0][0][0].alpha1rho1, Block_t::gptfloats); //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
                }
            }
        }
    };


    template < typename Kernel >
    struct Update_State
    {
        BlockInfo * ary;
        Real min_G, min_P, min_r;
        Real alpha, beta;

    public:
        Update_State(BlockInfo * ary, const Real min_G, const Real min_P, const Real min_r, const Real alpha, const Real beta): ary(ary), min_G(min_G), min_P(min_P), min_r(min_r), alpha(alpha), beta(beta) { }

        Update_State(const Update_State& c): ary(c.ary), min_G(c.min_G), min_P(c.min_P), min_r(c.min_r), alpha(c.alpha), beta(c.beta) { }

        void omp(const int N)
        {
#pragma omp parallel
            {
                /* TODO: (fabianw; Sun 25 Oct 2015 03:00:39 PM CET) testing
                 * this kernel */
                /* Kernel kernel(min_G,min_P, min_r, alpha, beta); */
                /* Kernel kernel(Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, Simulation_Environment::PC1, Simulation_Environment::PC2); */
                Kernel kernel(Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, Simulation_Environment::PC1, Simulation_Environment::PC2, state_alpha, state_beta);

#pragma omp for schedule(runtime)
                for(int r=0; r<N; ++r)
                {
                    Block_t & block = *(Block_t *)ary[r].ptrBlock;
                    kernel.compute(&block.data[0][0][0].alpha1rho1, Block_t::gptfloats); //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
                }
            }
        }
    };

    struct Update_Sponge
    {
        Real p_inf;
        BlockInfo * ary;
        FluidElement p;
        int s[3], e[3];

    public:

// TODO: VOLFRAC_EQSYS_URSULA: include again: old sponge
/*
        Update_Sponge(float p_inf, BlockInfo * ary): p_inf(p_inf), ary(ary)
        {
            const double pressure = p_inf;
            const double G1 = Simulation_Environment::GAMMA1-1;
            const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
            const double mix_gamma = 1 + G1;
            const double mix_pinf  = (mix_gamma-1)/mix_gamma * (F1/G1);
            p.G  = 1./(mix_gamma-1);
            p.P = mix_gamma*mix_pinf/(mix_gamma-1);
            p.energy   = pressure*p.G + p.P;
        }

        Update_Sponge(const Update_Sponge& c): p_inf(c.p_inf), ary(c.ary)
        {
            p.rho = 1000.0;
            const double pressure = p_inf;
            const double G1 = Simulation_Environment::GAMMA1-1;
            const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
            const double mix_gamma = 1 + G1;
            const double mix_pinf  = (mix_gamma-1)/mix_gamma * (F1/G1);
            p.G  = 1./(mix_gamma-1);
            p.P = mix_gamma*mix_pinf/(mix_gamma-1);
            p.energy   = pressure*p.G + p.P;
        }
*/

        Update_Sponge(float p_inf, BlockInfo * ary): p_inf(p_inf), ary(ary)
        {
            const double pressure = p_inf;
            const double G1 = Simulation_Environment::GAMMA1-1;
            const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
            p.energy = (pressure + F1) / G1;
        }

        Update_Sponge(const Update_Sponge& c): p_inf(c.p_inf), ary(c.ary)
        {
            const double pressure = p_inf;
            const double G1 = Simulation_Environment::GAMMA1-1;
            const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
            p.energy = (pressure + F1) / G1;
        }

        /* TODO: (fabianw; Tue 29 Sep 2015 06:50:27 PM CEST) these are the OG
         * two methods */
        /* template<int dir, int side> */
        /* void _setup_sponge() */
        /* { */
        /*     s[0] =	s[1] =	s[2] =	0; */
        /*     e[0] =  Block_t::sizeX; */
        /*     e[1] =  Block_t::sizeY; */
        /*     e[2] =	Block_t::sizeZ; */
        /* } */

        /* template<int dir, int side> */
        /* void applyBC_sponge_subsonic(Block_t & block, const Real h) */
        /* { */
        /*     _setup_sponge<dir,side>(); */

        /*     for(int iz=s[2]; iz<e[2]; iz++) */
        /*         for(int iy=s[1]; iy<e[1]; iy++) */
        /*             for(int ix=s[0]; ix<e[0]; ix++) */
        /*             { */
        /*                 //extend the last or first point all over the first or last block */
        /*                 block(ix,iy,iz) = block(dir==0? (side==0? Block_t::sizeX-1:0):ix, */
        /*                                         dir==1? (side==0? Block_t::sizeY-1:0):iy, */
        /*                                         dir==2? (side==0? Block_t::sizeZ-1:0):iz); */

        /*                 const Real ke = 0.5*(pow(block(ix,iy,iz).u,2)+pow(block(ix,iy,iz).v,2)+pow(block(ix,iy,iz).w,2))/block(ix,iy,iz).rho; */

        /*                 const int offset = dir==0? (side==0? Block_t::sizeX-1-ix : ix) : dir==1? (side==0? Block_t::sizeY-1-iy : iy) : (side==0? Block_t::sizeZ-1-iz : iz); */

        /*                 const int extent = dir==0? Block_t::sizeX : dir==1? Block_t::sizeY : Block_t::sizeZ; */

        /*                 const Real EPSILON = 0.5*(double)extent; */
        /*                 const Real phi = (double)offset; */
        /*                 const Real alpha = M_PI*min(1., max(0., 0.5*phi/EPSILON)); */
        /*                 const Real blender = 0.5-0.5*cos(alpha); */

        /*                 const Real difference = p.energy+ke-block(dir==0? (side==0? Block_t::sizeX-1:0):ix, */
        /*                                                        dir==1? (side==0? Block_t::sizeY-1:0):iy, */
        /*                                                        dir==2? (side==0? Block_t::sizeZ-1:0):iz).energy; */

        /*                 block(ix,iy,iz).energy = difference*(double)offset*blender/(double)(extent) + block(dir==0? (side==0? Block_t::sizeX-1:0):ix, */
        /*                                                                                                     dir==1? (side==0? Block_t::sizeY-1:0):iy, */
        /*                                                                                                     dir==2? (side==0? Block_t::sizeZ-1:0):iz).energy; */
        /*             } */

        /* }; */
        /* TODO: (fabianw; Tue 29 Sep 2015 06:50:44 PM CEST) end og sponge
         * methods */

        template<int dir, int side>
        void _setup_sponge()
        {
            const int width_x = min (LSRK3data::sponge_width [0], Block_t::sizeX);
            const int width_y = min (LSRK3data::sponge_width [1], Block_t::sizeY);
            const int width_z = min (LSRK3data::sponge_width [2], Block_t::sizeZ);

            s[0] = dir == 0 ? (side == 0 ? 0 : Block_t::sizeX - width_x) : 0;
            s[1] = dir == 1 ? (side == 0 ? 0 : Block_t::sizeY - width_y) : 0;
            s[2] = dir == 2 ? (side == 0 ? 0 : Block_t::sizeZ - width_z) : 0;

            e[0] = dir == 0 ? (side == 0 ? width_x : Block_t::sizeX) : Block_t::sizeX;
            e[1] = dir == 1 ? (side == 0 ? width_y : Block_t::sizeY) : Block_t::sizeY;
            e[2] = dir == 2 ? (side == 0 ? width_z : Block_t::sizeZ) : Block_t::sizeZ;
        }

        template<int dir, int side>
        void applyBC_sponge_subsonic(Block_t & block, const Real h)
        {
            _setup_sponge<dir,side>();

            for(int iz=s[2]; iz<e[2]; iz++)
                for(int iy=s[1]; iy<e[1]; iy++)
                    for(int ix=s[0]; ix<e[0]; ix++)
                    {
                        // extend the last or first cell in the specified direction over the specified extent
                        block (ix, iy, iz) = block (dir==0 ? (side==0 ? e[0]-1 : s[0]) : ix,
                                                    dir==1 ? (side==0 ? e[1]-1 : s[1]) : iy,
                                                    dir==2 ? (side==0 ? e[2]-1 : s[2]) : iz);

                        //const Real ke = 0.5*(pow(block(ix,iy,iz).u,2)+pow(block(ix,iy,iz).v,2)+pow(block(ix,iy,iz).w,2))/block(ix,iy,iz).rho;

                        const int offset = dir==0 ? (side==0 ? e[0]-1-ix : ix-s[0]) : dir==1 ? (side==0? e[1]-1-iy : iy-s[1]) : (side==0 ? e[2]-1-iz : iz-s[2]);

                        const int extent = e[dir] - s[dir];

                        const Real EPSILON = 0.5 * (double)extent;
                        const Real phi = (double)offset;
                        const Real alpha = M_PI * min (1., max(0., 0.5*phi/EPSILON));
                        const Real blender = 0.5 - 0.5 * cos (alpha);

                        const Real energy = block (dir==0 ? (side==0 ? e[0]-1 : s[0]) : ix,
                                                   dir==1 ? (side==0 ? e[1]-1 : s[1]) : iy,
                                                   dir==2 ? (side==0 ? e[2]-1 : s[2]) : iz
                                                  ).energy;

                        const Real difference = p.energy - energy;

                        block(ix,iy,iz).energy += difference * (double)offset * blender / (double)extent;
                        //block(ix,iy,iz).energy = difference*(double)offset*blender/(double)(extent) + block(dir==0? (side==0? e[0]-1:s[0]):ix,
//                dir==1? (side==0? e[1]-1:s[1]):iy,
//                dir==2? (side==0? e[2]-1:s[2]):iz).energy;
                    }

        };

        void omp(const int N, const int NX, const int NY, const int NZ)
        {
#pragma omp parallel
            {
#ifdef _USE_NUMA_
                const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
                const int mynode = omp_get_thread_num() / cores_per_node;
                numa_run_on_node(mynode);
#endif

#pragma omp for schedule(runtime)
                for(int r=0; r<N; ++r)
                {
                    Block_t & block = *(Block_t *)ary[r].ptrBlock;
                    const double spacing = ary[r].h_gridpoint;

                    //sponge magic
#ifndef _BCLABCLOUDSYMABSORB_
                    if (ary[r].index[0]==0)       applyBC_sponge_subsonic<0,0>(block, spacing);
                    if (ary[r].index[0]==NX-1)    applyBC_sponge_subsonic<0,1>(block, spacing);
                    if (ary[r].index[1]==0)       applyBC_sponge_subsonic<1,0>(block, spacing);
                    if (ary[r].index[1]==NY-1)    applyBC_sponge_subsonic<1,1>(block, spacing);
                    if (ary[r].index[2]==0)       applyBC_sponge_subsonic<2,0>(block, spacing);
                    if (ary[r].index[2]==NZ-1)    applyBC_sponge_subsonic<2,1>(block, spacing);
#else
                    if (ary[r].index[0]==NX-1)    applyBC_sponge_subsonic<0,1>(block, spacing);
                    if (ary[r].index[1]==NY-1)    applyBC_sponge_subsonic<1,1>(block, spacing);
                    if (ary[r].index[2]==NZ-1)    applyBC_sponge_subsonic<2,1>(block, spacing);
#endif
                    /*
                    if (ary[r].index[0]==0 && ! LSRK3data::BC_PERIODIC[0])             applyBC_sponge_subsonic<0,0>(block, spacing);
                    if (ary[r].index[0]==NX-1 && ! LSRK3data::BC_PERIODIC[0])          applyBC_sponge_subsonic<0,1>(block, spacing);
                    if (ary[r].index[1]==0 && ! LSRK3data::BC_PERIODIC[1])             applyBC_sponge_subsonic<1,0>(block, spacing);
                    if (ary[r].index[1]==NY-1 && ! LSRK3data::BC_PERIODIC[1])    applyBC_sponge_subsonic<1,1>(block, spacing);
                    if (ary[r].index[2]==0 && ! LSRK3data::BC_PERIODIC[2])             applyBC_sponge_subsonic<2,0>(block, spacing);
                    if (ary[r].index[2]==NZ-1 && ! LSRK3data::BC_PERIODIC[2])    applyBC_sponge_subsonic<2,1>(block, spacing);
                    */
                }
            }
        }
    };
};


class FlowStep_LSRK3
{
private:
    Grid_t& grid;

protected:
    Real smoothlength, h_min;
    Real PEAKPERF_CORE, PEAKBAND;
    string blockdispatcher;
    const int verbosity;

    Real _computeSOS();

    ArgumentParser& parser;

public:
    Real CFL;

    FlowStep_LSRK3(Grid_t& grid, const Real CFL, ArgumentParser& parser, const int verbosity=1):
    grid(grid), CFL(CFL), parser(parser), verbosity(verbosity)
    {
        parser.unset_strict_mode();

        PEAKPERF_CORE = parser("-pp").asDouble(27.2*2);
        PEAKBAND = parser("-pb").asDouble(19);
        blockdispatcher = parser("-dispatcher").asString("none");

        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        h_min = vInfo[0].h_gridpoint;
        // TODO: [fabianw@mavt.ethz.ch; Wed May 10 2017 03:56:06 PM (-0700)]
        // this must be localized to grid spacing
        smoothlength = (Real)(parser("-mollfactor").asInt())*sqrt(3.)*h_min;
        Simulation_Environment::EPSILON = smoothlength;
    }

    Real operator()(const Real max_dt, const Real current_time=0.0);

    inline void set_CFL(const Real _CFL) {CFL = _CFL;}
    inline void set_hmin(const Real _h)
    {
        h_min = _h;
        smoothlength = (Real)(parser("-mollfactor").asInt())*sqrt(3.)*h_min;
        Simulation_Environment::EPSILON = smoothlength;
    }
    inline Real get_hmin() const { return h_min; }
    void set_constants();
};
