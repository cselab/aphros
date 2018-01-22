/*
 *  FlowStep_LSRK3.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include "FlowStep_LSRK3.h"

#include <cstdlib>
#include <omp.h>

using namespace std;

namespace LSRK3data
{
    Real smoothlength = -1;
    int verbosity;
    float PEAKPERF_CORE, PEAKBAND;
    int  NCORES;
    string dispatcher;
    int ReportFreq = 1;
    int step_id = 0;

    int sponge = 0;
    Real sponge_pref;
    int sponge_width [3] = {0, 0, 0};

    int state = 0;
    Real state_alpha = 2.0;
    Real state_beta  = 4.0;
}


template < typename TSOS>
Real _computeSOS_OMP(Grid_t& grid)
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

    Real * const local_sos = tmp;

#pragma omp parallel
    {

#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif

        TSOS kernel(Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2,
                    Simulation_Environment::PC1, Simulation_Environment::PC2);
#pragma omp for schedule(runtime)
        for (int i=0; i<N; ++i)
        {
            Block_t & block = *(Block_t *)ary[i].ptrBlock;
            local_sos[i] =  kernel.compute(&block.data[0][0][0].alpha1rho1, Block_t::gptfloats); //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
        }
    }

    Real global_sos = local_sos[0];

#pragma omp parallel
    {
        Real mymax = local_sos[0];

#pragma omp for schedule(runtime)
        for(int i=0; i<N; ++i)
            mymax = max(local_sos[i], mymax);

#pragma omp critical
        {
            global_sos = max(global_sos, mymax);
        }
    }

    free(tmp);

#else   /* FUSED_MAXSOS */

    Real global_sos = 0;

#pragma omp parallel
    {
        /* TSOS kernel; */
        TSOS kernel(Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2,
                    Simulation_Environment::PC1, Simulation_Environment::PC2);

#pragma omp for schedule(runtime) reduction(max:global_sos)
        for (size_t i=0; i<N; ++i)
        {
            Block_t & block = *(Block_t *)ary[i].ptrBlock;
            global_sos =  max(global_sos, kernel.compute(&block.data[0][0][0].alpha1rho1, Block_t::gptfloats)); //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
        }
    }

#endif

#ifdef _USE_HPM_
    HPM_Stop("dt");
#endif

    return global_sos;
}


Real FlowStep_LSRK3::_computeSOS()
{
    Real sos = -1;

    const string kernels = parser("-kernels").asString("cpp");
    vector<BlockInfo> vInfo = grid.getBlocksInfo();

    Timer timer;

    timer.start();

#if defined(_QPX_) || defined(_QPXEMU_)
    if (kernels == "qpx")
        sos = _computeSOS_OMP<MaxSpeedOfSound_QPX_5eq>(grid);
    else
#endif
        sos = _computeSOS_OMP<MaxSpeedOfSound_CPP_5eq>(grid);

    const Real time = timer.stop();

    if (LSRK3data::verbosity >= 1 && LSRK3data::step_id % LSRK3data::ReportFreq == 0)
    {
        MaxSpeedOfSound_CPP::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), time);

        cout << "MAXSOS: " << time << "s (per substep), " << time/vInfo.size()*1e3 << " ms (per block)" << endl;
    }

    return sos;
}


template<typename Kflow, typename Kdiff, typename Kupdate, typename Ksharp, typename Ksource>
struct LSRKstep
{
    template<typename T>
    void check_error(const double tol, T ref[], T val[], const int N)
    {
        for(int i=0; i<N; ++i)
        {
            assert(!std::isnan(ref[i]));

            assert(!std::isnan(val[i]));

            const double err = ref[i] - val[i];

            const double relerr = err/std::max(1e-6, (double)std::max(fabs(val[i]), fabs(ref[i])));

            if (LSRK3data::verbosity>= 1) printf("+%1.1e,", relerr);

            if (fabs(relerr) >= tol && fabs(err) >= tol)
                printf("\n%d: %e %e -> %e %e\n", i, ref[i], val[i], err, relerr);

            //assert(fabs(relerr) < tol || fabs(err) < tol);
        }
        if (LSRK3data::verbosity >=1) printf("\t");
    }

    template<typename T>
    void check_error(const double tol, T ref, T val)
    {
        check_error(tol, &ref, &val, 1);
    }

    template<typename T>
    void check_error(const double tol, T ref[], T val[], const int N, const int bidx[], const int idx[], const int comp)
    {
        for(int i=0; i<N; ++i)
        {
            assert(!std::isnan(ref[i]));

            assert(!std::isnan(val[i]));

            const double err = fabs(ref[i]) - fabs(val[i]);

            const double relerr = err/std::max(1e-6, (double)std::max(fabs(val[i]), fabs(ref[i])));

            if (LSRK3data::verbosity>= 1) printf("+%1.1e,", relerr);

            if (fabs(relerr) >= tol && fabs(err) >= tol)
                printf("\n%d: %e %e -> %e %e, located at mirror block %d,%d,%d and point %d,%d,%d, comp %d \n", i, ref[i], val[i], err, relerr, bidx[0], bidx[1], bidx[2], idx[0], idx[1], idx[2], comp);


            //assert(fabs(relerr) < tol || fabs(err) < tol);
        }
        if (LSRK3data::verbosity >=1) printf("\t");
    }

    // TODO: URSULA this function seems to be unused
    void _check_symmetry(Grid_t& grid)
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();

        for(int i=0; i<(int)vInfo.size(); i++)
        {
            BlockInfo info = vInfo[i];

            if (info.index[1]>grid.getBlocksPerDimension(1)/2-1) continue;

            Block_t& b = *(Block_t*)info.ptrBlock;

            const int bidx_mirror[3] = {info.index[0], grid.getBlocksPerDimension(1)-info.index[1]-1, info.index[2]};
            const int i_mirror = bidx_mirror[0] + (bidx_mirror[1] + bidx_mirror[2]*grid.getBlocksPerDimension(1))*grid.getBlocksPerDimension(0);
            assert(i_mirror<grid.getBlocksPerDimension(0)*grid.getBlocksPerDimension(1)*grid.getBlocksPerDimension(2) && i_mirror>=0);

            BlockInfo info_mirror = vInfo[i_mirror];
            Block_t& b_mirror = *(Block_t*)info_mirror.ptrBlock;

            for(int iz=0; iz<Block_t::sizeZ; iz++)
                for(int iy=0; iy<Block_t::sizeY; iy++)
                    for(int ix=0; ix<Block_t::sizeX; ix++)
                    {
                        const int idx_mirror[3] = {ix, Block_t::sizeY-iy-1, iz};

                        check_error(std::numeric_limits<Real>::epsilon()*50, &b(ix,iy,iz).alpha1rho1, &b_mirror(idx_mirror[0],idx_mirror[1],idx_mirror[2]).alpha1rho1, 7, bidx_mirror, idx_mirror, 0); //TODO: URSULA: This may easily be generalized with an option to address the first element. Name does not matter here.
                        //                        check_error(std::numeric_limits<Real>::epsilon()*50, &b.tmp[iz][iy][ix][0], &b_mirror.tmp[idx_mirror[2]][idx_mirror[1]][idx_mirror[0]][0], 7, bidx_mirror, idx_mirror, 1);
                    }
        }

        cout << "Symmetry check done" << endl;
    }

    LSRKstep(Grid_t& grid, const Real dt, const Real current_time, const double epsilonSharpening=-1)
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();

        vector< vector<double> > timings;

        //_check_symmetry(grid);

        double A[3], B[3];

        // Williamson
        A[0] = 0.0;     B[0] = 1./4;
        A[1] = -17./32; B[1] = 8./9;
        A[2] = -32./27; B[2] = 3./4;

        // Gottlieg & Shu
        // A[0] = 0.0;                B[0] = 0.924574000000000;
        // A[1] = -2.915492524638791; B[1] = 0.287713063186749;
        // A[2] = -0.000000093517376; B[2] = 0.626538109512740;

        timings.push_back(step(grid, vInfo, A[0], B[0], dt, current_time, epsilonSharpening));
        timings.push_back(step(grid, vInfo, A[1], B[1], dt, current_time + B[0]*dt, epsilonSharpening));
        timings.push_back(step(grid, vInfo, A[2], B[2], dt, current_time + (B[0]+(A[1]+1.0)*B[1])*dt, epsilonSharpening));

        // forward Euler
        // timings.push_back(step(grid, vInfo, 0.0, 1.0, dt, current_time, epsilonSharpening));

        const double avg1 = ( timings[0][0] + timings[1][0] + timings[2][0] )/3;
        const double avg2 = ( timings[0][1] + timings[1][1] + timings[2][1] )/3;

        if (LSRK3data::verbosity >= 1 && LSRK3data::step_id % LSRK3data::ReportFreq == 0)
        {
            cout << "FLOWSTEP: " << avg1 << "s (per substep), " << avg1/vInfo.size()*1e3 << " ms (per block)" << endl;
            cout << "UPDATE: " << avg2 << "s (per substep), " << avg2/vInfo.size()*1e3 << " ms (per block)" << endl;
            //const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NTIMES, const int NBLOCKS, const float MEASUREDTIME
            Kflow::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg1);
            Kupdate::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg2);

           // TODO: Ursula (October 17 2015): older vesions of Cubism-MPCF considered also diffusion and surface tension
        }
    }

    vector<double> step(Grid_t& grid, vector<BlockInfo>& vInfo, Real a, Real b, const Real dt, const Real current_time, const double epsilonSharpening)
    {
        Timer timer;
        vector<double> res;

#ifndef _NONUNIFORM_BLOCK_
        const Real h = vInfo[0].h_gridpoint; // uniform grid size
#endif /* _NONUNIFORM_BLOCK_ */

///////////////////////////////////////////////////////////////////////////////
// Stencil (Lab) based schemes
///////////////////////////////////////////////////////////////////////////////
#ifdef _USE_HPM_
        HPM_Start("RHS");
#endif

#ifdef _NONUNIFORM_BLOCK_
        // nonuniform
        LSRK3data::FlowStepNonuniform<Kflow, Lab> rhs(a, dt);
        timer.start();
        processOMP<Lab>(rhs, grid, current_time);
        const double t1 = timer.stop();
#else
        // uniform
        LSRK3data::FlowStep<Kflow, Lab> rhs(a, dt/h);
        timer.start();
        processOMP<Lab>(rhs, grid, current_time);
        const double t1 = timer.stop();
#endif /* _NONUNIFORM_BLOCK_ */

#ifdef _USE_HPM_
        HPM_Stop("RHS");
#endif
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#ifdef _USE_HPM_
    if (LSRK3data::step_id>0)             HPM_Start("DIFFUSION");
#endif
    if(Simulation_Environment::MU_MAX>0)
    {
#ifdef _NONUNIFORM_BLOCK_
        // TODO: [fabianw@mavt.ethz.ch; Wed May 10 2017 07:11:27 PM (-0700)]
#else
        LSRK3data::Diffusion<Kdiff, Lab> diffusion(dt/h, Simulation_Environment::MU1, Simulation_Environment::MU2);
        timer.start();
        processOMP<Lab>(diffusion, grid, current_time);
        const double tdiffusion = timer.stop();
#endif /* _NONUNIFORM_BLOCK_ */
    }
#ifdef _USE_HPM_
    if (LSRK3data::step_id>0)             HPM_Stop("DIFFUSION");
#endif
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#ifdef _USE_HPM_
    HPM_Start("SHARPENING");
#endif

    if(epsilonSharpening>0.0)
    {
#ifdef _NONUNIFORM_BLOCK_
        // TODO: [fabianw@mavt.ethz.ch; Wed May 10 2017 07:13:33 PM (-0700)]
#else
        const Real U0 = computeMaxIVel(grid);
        LSRK3data::InterfaceSharpening<Ksharp, Lab> intersharp(dt/h, epsilonSharpening, U0);
        processOMP<Lab>(intersharp, grid, current_time);
#endif /* _NONUNIFORM_BLOCK_ */
    }

#ifdef _USE_HPM_
    HPM_Stop("SHARPENING");
#endif
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Pointwise schemes
///////////////////////////////////////////////////////////////////////////////
#ifdef _USE_HPM_
    if (LSRK3data::step_id>0)             HPM_Start("SOURCE");
#endif
    // FIXME: [fabianw@mavt.ethz.ch; Thu May 11 2017 09:41:27 PM (-0700)]
    const bool source = true; // for now

    if(source)
    {
        // TODO: [fabianw@mavt.ethz.ch; Wed May 17 2017 03:47:58 PM (-0700)]
        // parameter need to be specified in a more convenient way
        OneWayAcousticSource_CPP::SourceParameter p;
        p.gamma = Simulation_Environment::GAMMA1;
        p.c0 = Simulation_Environment::C1;
        p.x0 = 0.5*Simulation_Environment::extents[0];
        p.smooth = Simulation_Environment::EPSILON;

        p.sigma_t = 5.0e-6;
        p.t0 = 20.0e-6;
        p.amplitude = 100*1.0e-5;
        // p.sigma_t = 1.5;
        // p.t0 = 4.5;
        // p.amplitude = 1.0e-5;
        // p.amplitude = 1.0;

        LSRK3data::OneWayAcousticSource<OneWayAcousticSource_CPP> source(current_time, dt, p,  &vInfo.front());

        // LSRK3data::Source<Ksource> source (1.0, dt, Simulation_Environment::grav_accel, Simulation_Environment::vol_bodyforce,  &vInfo.front());

        source.omp(vInfo.size());
    }
#ifdef _USE_HPM_
    if (LSRK3data::step_id>0)            HPM_Stop("SOURCE");
#endif
///////////////////////////////////////////////////////////////////////////////


    if (Simulation_Environment::SIGMA > 0.0)
        cout << "WARNING: SURFACE TESION NOT YET INCLUDED IN TIME STEPPING OF NODE LEVEL!" << endl;


    // TODO: [fabianw@mavt.ethz.ch; Wed May 10 2017 07:16:41 PM (-0700)] does
    // this need to be here? Need to check this
#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
    //1d characteristic-based bc currently only for Euler equations
    LSRK3data::UpdateBC<Lab> upbc(dt,a,b);
    processOMP<Lab>(upbc, grid, current_time);
#endif

///////////////////////////////////////////////////////////////////////////////
// RK update
///////////////////////////////////////////////////////////////////////////////
#ifdef _USE_HPM_
        HPM_Start("Update");
#endif
        LSRK3data::Update<Kupdate> update(b, &vInfo.front());

        timer.start();
        update.omp(vInfo.size());
        const double t2 = timer.stop();
#ifdef _USE_HPM_
        HPM_Stop("Update");
#endif
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#ifdef _USE_HPM_
    if (LSRK3data::step_id>0)             HPM_Start("Update State");
#endif
    if(LSRK3data::state==1)
    {
        // std::cout << "######################################\n#      UPDATE STATE DEACTIVATED      #\n######################################" << std::endl;
        // LSRK3data::Update_State<Update_State_CPP_5eq> update_state(&vInfo.front(), Simulation_Environment::min_G, Simulation_Environment::min_P, 0.0, LSRK3data::state_alpha, LSRK3data::state_beta);
        // update_state.omp(vInfo.size());
    }
#ifdef _USE_HPM_
    if (LSRK3data::step_id>0)             HPM_Stop("Update State");
#endif
///////////////////////////////////////////////////////////////////////////////

        res.push_back(t1);
        res.push_back(t2);

        return res;
    }

    Real computeMaxIVel(Grid_t& grid)
    {
        Real ivel = -1;

        Timer timer;

        timer.start();
        ivel = _computeMaxIVel_OMP<MaxInterfaceVel_CPP_5eq>(grid);
        const Real time = timer.stop();

        return ivel;
    }
};

void FlowStep_LSRK3::set_constants()
{
    LSRK3data::smoothlength = smoothlength;
    LSRK3data::verbosity = verbosity;
    LSRK3data::PEAKBAND = PEAKBAND;
    LSRK3data::PEAKPERF_CORE = PEAKPERF_CORE;
    LSRK3data::NCORES = parser("-ncores").asInt(1);
    LSRK3data::dispatcher = blockdispatcher;
    LSRK3data::ReportFreq = parser("-report").asInt(50);
    LSRK3data::sponge = parser("-sponge").asInt(0);
    LSRK3data::sponge_pref = Simulation_Environment::P1;
    LSRK3data::sponge_width [0] = parser("-spongewidth").asInt(32);
    LSRK3data::sponge_width [1] = parser("-spongewidth").asInt(32);
    LSRK3data::sponge_width [2] = parser("-spongewidth").asInt(32);
    LSRK3data::state = parser("-state").asInt(0);
    LSRK3data::state_alpha = parser("-state_alpha").asDouble(2.0);
    LSRK3data::state_beta  = parser("-state_beta").asDouble(4.0);
}

Real FlowStep_LSRK3::operator()(const Real max_dt, const Real current_time)
{
    set_constants();

    Timer timer;
    timer.start();
#ifdef _USE_HPM_
    HPM_Start("dt");
#endif

    const Real maxSOS = _computeSOS();

#ifdef _USE_HPM_
    HPM_Stop("dt");
#endif

    const double t_sos = timer.stop();

    cout << "sos take " << t_sos << " sec" << endl;

    double dt = min(max_dt, CFL*h_min/maxSOS);
    cout << "sos max is " << setprecision(8) << maxSOS << ", " << "dt is "<< dt << "\n";

    if (maxSOS > 1e6)
    {
        cout << "Speed of sound is too high. Is it realistic?" << endl;
        abort();
    }

    if (dt<std::numeric_limits<double>::epsilon()*1e1)
    {
        cout << "Last time step encountered." << endl;
        return 0;
    }

    if (LSRK3data::verbosity >= 1)
        cout << "Dispatcher is " << LSRK3data::dispatcher << endl;


    // process step
    const double epsilonSharpening = parser("-eps-sharp").asDouble(-1.0);
    const bool compact = parser("-sharp-form").asBool(true);

#ifdef _NONUNIFORM_BLOCK_
    if (parser("-kernels").asString("cpp")=="cpp")
        LSRKstep<Convection_CPP_HLLC_5eq_nonuniform, Diffusion_CPP, Update_CPP, InterfaceSharpening_CPP_cell_compact_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
// #if defined(_QPX_) || defined(_QPXEMU_)
//     else if (parser("-kernels").asString("cpp")=="qpx")
// #if defined(_BGQ_)
//         LSRKstep<Convection_QPX_HLLC_5eq, Diffusion_CPP, Update_CPP, InterfaceSharpening_CPP_cell_compact_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
// #else
//     LSRKstep<Convection_QPX_HLLC_5eq, Diffusion_CPP, Update_QPX, InterfaceSharpening_CPP_cell_compact_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
// #endif
// #endif
    else
    {
        cout << "combination not supported yet" << endl;
        abort();
    }

#else

    // TODO: [fabianw@mavt.ethz.ch; Wed May 10 2017 03:42:49 PM (-0700)] needs
    // cleaning...
    //
    if (compact)
    {
        if (parser("-kernels").asString("cpp")=="cpp")
            LSRKstep<Convection_CPP_HLLC_5eq, Diffusion_CPP, Update_CPP, InterfaceSharpening_CPP_cell_compact_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
#if defined(_QPX_) || defined(_QPXEMU_)
        else if (parser("-kernels").asString("cpp")=="qpx")
#if defined(_BGQ_)
            LSRKstep<Convection_QPX_HLLC_5eq, Diffusion_CPP, Update_CPP, InterfaceSharpening_CPP_cell_compact_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
#else
            LSRKstep<Convection_QPX_HLLC_5eq, Diffusion_CPP, Update_QPX, InterfaceSharpening_CPP_cell_compact_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
#endif
#endif
        else
        {
            cout << "combination not supported yet" << endl;
            abort();
        }
    }
    else
    {
        if (parser("-kernels").asString("cpp")=="cpp")
            LSRKstep<Convection_CPP_HLLC_5eq, Diffusion_CPP, Update_CPP, InterfaceSharpening_CPP_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
#if defined(_QPX_) || defined(_QPXEMU_)
        else if (parser("-kernels").asString("cpp")=="qpx")
#if defined(_BGQ_)
            LSRKstep<Convection_QPX_HLLC_5eq, Diffusion_CPP, Update_CPP, InterfaceSharpening_CPP_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
#else
            LSRKstep<Convection_QPX_HLLC_5eq, Diffusion_CPP, Update_QPX, InterfaceSharpening_CPP_5eq, Source_CPP>(grid, dt, current_time, epsilonSharpening);
#endif
#endif
        else
        {
            cout << "combination not supported yet" << endl;
            abort();
        }
    }

#endif /* _NONUNIFORM_BLOCK_ */

    LSRK3data::step_id++;

    return dt;
}


