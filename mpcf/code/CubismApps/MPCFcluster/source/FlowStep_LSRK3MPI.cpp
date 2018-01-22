/*
 *  FlowStep_LSRK3MPI.cpp
 *  MPCFcluster
 *
 *  Created by Babak Hejazi on 11/20/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <mpi.h>
#include "FlowStep_LSRK3MPI.h"

using namespace std;

namespace LSRK3MPIdata
{
    double t_fs = 0;
    double t_up = 0;
    double t_synch_fs = 0;
    double t_bp_fs = 0;
    int counter = 0;
    int GSYNCH = 0;
    int nsynch = 0;

#ifndef _SEQUOIA_
    MPI_ParIO_Group hist_group;
    MPI_ParIO hist_update, hist_rhs, hist_stepid, hist_nsync;
#endif

    template<typename Kflow, typename Kupdate>
    void notify(double avg_time_rhs, double avg_time_update, const size_t NBLOCKS, const size_t NTIMES, MPI_Comm comm)
    {
#ifndef _SEQUOIA_
        if(LSRK3data::step_id % LSRK3data::ReportFreq == 0 && LSRK3data::step_id > 0)
        {
            hist_update.Consolidate(LSRK3data::step_id);
            hist_stepid.Consolidate(LSRK3data::step_id);
            hist_rhs.Consolidate(LSRK3data::step_id);
            hist_nsync.Consolidate(LSRK3data::step_id);
        }

        hist_update.Notify((float)avg_time_update);
        hist_rhs.Notify((float)avg_time_rhs);
        hist_stepid.Notify((float)LSRK3data::step_id);
        hist_nsync.Notify((float)nsynch/NTIMES);
#endif

        nsynch = 0;

        if(LSRK3data::step_id % LSRK3data::ReportFreq == 0 && LSRK3data::step_id > 0)
        {
            double global_t_synch_fs = 0, global_t_bp_fs = 0;
            double global_avg_time_rhs = 0, global_avg_time_update = 0;
            double global_t_fs, global_t_up;

            int global_counter = 0;

            MPI_Reduce(&t_synch_fs, &global_t_synch_fs, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&t_bp_fs, &global_t_bp_fs, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&counter, &global_counter, 1, MPI_INT, MPI_SUM, 0, comm);
            MPI_Reduce(&avg_time_rhs, &global_avg_time_rhs, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&avg_time_update, &global_avg_time_update, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&t_fs, &global_t_fs, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&t_up, &global_t_up, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

            t_synch_fs = t_bp_fs = t_fs = t_up = counter = 0;

            global_t_synch_fs /= NTIMES;
            global_t_bp_fs /= NTIMES;
            global_counter /= NTIMES;
            global_t_fs /= NTIMES;
            global_t_up /= NTIMES;

            int NRANKS;
            MPI_Comm_size(comm, &NRANKS);

            if (LSRK3data::verbosity >= 1)
            {
                cout << "FLOWSTEP: " << avg_time_rhs << " s (per substep), " << avg_time_rhs/NBLOCKS*1e3 << " ms (per block) " << global_avg_time_rhs/NBLOCKS*1e3/NRANKS << " ms (per block per node)" << endl;

                cout << "TIME LOCALLY AVERAGED FLOWSTEP: " << global_avg_time_rhs/NRANKS <<" s (per substep per node), " << global_avg_time_rhs/NBLOCKS*1e3/NRANKS << " ms (per substep per node per block)" << endl;

                cout << "TIME GLOBALLY AVERAGED FLOWSTEP: " << global_t_fs/NRANKS/(double)LSRK3data::ReportFreq << " s (per substep)" << endl;

                cout << "===========================STAGE===========================" << endl;
                cout << "Synch done in "<< global_counter/NRANKS/(double)LSRK3data::ReportFreq << " passes" << endl;
                cout << "SYNCHRONIZER FLOWSTEP "<< global_t_synch_fs/NRANKS/(double)LSRK3data::ReportFreq << " s" << endl;
                cout << "BP FLOWSTEP "<< global_t_bp_fs/NRANKS/(double)LSRK3data::ReportFreq << " s" << endl;
                cout << "======================================================" << endl;

                Kflow::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1,  NBLOCKS*NRANKS, global_t_fs/(double)LSRK3data::ReportFreq/NRANKS);
                Kupdate::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, NBLOCKS*NRANKS, global_t_up/(double)LSRK3data::ReportFreq/NRANKS);

               // TODO: Ursula (October 17 2015) old cubism version included viscous and surface tension term here as well

            }
        }
    }
}

template<typename TGrid>
template<typename Kflow, typename Ksource, typename Kdiff, typename Ksurf, typename Kupdate, typename Kupdate_state, typename Ksharp>
pair<Real, Real> FlowStep_LSRK3MPI<TGrid>::step(TGrid& grid, vector<BlockInfo>& vInfo, const Real a, const Real b, const Real dt, const Real current_time)
{
    const bool record = LSRK3data::step_id%10==0;

    const Real h = vInfo[0].h_gridpoint; // uniform grid size
    const Real dtinvh = dt / h;

    // diffusion
    LSRK3data::Diffusion<Kdiff, Lab> diffusion(dtinvh, Simulation_Environment::MU1, Simulation_Environment::MU2);
    process< LabMPI >(diffusion, (TGrid&)grid, current_time, record);

    // update
    LSRK3data::Update<Kupdate> update(b, &vInfo.front());
    update.omp(vInfo.size());

    return pair<double, double>(1.,1.);
}

template<typename TGrid>
template<typename Kflow, typename Ksource, typename Kdiff, typename Ksurf, typename Kupdate, typename Kupdate_state, typename Ksharp>
inline void FlowStep_LSRK3MPI<TGrid>::_process_LSRK3(TGrid& grid, const Real dt, const Real current_time)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();

    vector< pair<double, double> > timings;

    string LSRK3coeffs = parser("-LSRK3coeffs").asString("williamson");

    double A [3], B [3];

    if (LSRK3coeffs == "williamson") {
      A [0] = 0.0;     B [0] = 1./4;
      A [1] = -17./32; B [1] = 8./9;
      A [2] = -32./27; B [2] = 3./4;
    }

    else if (LSRK3coeffs == "gottlieb-shu") {
      A [0] = 0.0;                B [0] = 0.924574000000000;
      A [1] = -2.915492524638791; B [1] = 0.287713063186749;
      A [2] = -0.000000093517376; B [2] = 0.626538109512740;
    }

    else {
      cout << " :: ERROR: Unknown LSRK3 coefficiets \"" << LSRK3coeffs << "\"" << endl;
      exit(-1);
    }

    timings.push_back(step<Kflow, Ksource, Kdiff, Ksurf, Kupdate, Kupdate_state, Ksharp>(grid, vInfo, A[0], B[0], dt, current_time));
    timings.push_back(step<Kflow, Ksource, Kdiff, Ksurf, Kupdate, Kupdate_state, Ksharp>(grid, vInfo, A[1], B[1], dt, current_time + B[0]*dt));
    timings.push_back(step<Kflow, Ksource, Kdiff, Ksurf, Kupdate, Kupdate_state, Ksharp>(grid, vInfo, A[2], B[2], dt, current_time + (B[0]+(A[1]+1.0)*B[1])*dt));

    double avg1 = ( timings[0].first  + timings[1].first  + timings[2].first  )/3;
    double avg2 = ( timings[0].second + timings[1].second + timings[2].second )/3;

    LSRK3MPIdata::notify<Kflow, Kupdate>(avg1, avg2, vInfo.size(), 3, grid.getCartComm());
}

template<typename TGrid>
Real FlowStep_LSRK3MPI<TGrid>::_computeSOS()
{
    double maxSOS;

    double local_maxSOS = FlowStep_LSRK3::_computeSOS();

    MPI_Comm mycart = grid.getCartComm();

    MPI_Allreduce(&local_maxSOS, &maxSOS, 1, MPI_DOUBLE, MPI_MAX, mycart);

    return maxSOS;
}

// // for interface sharpening: compute maximal interface velocity
// // TODO: (fabianw@mavt.ethz.ch; Wed 06 Apr 2016 11:26:22 AM CEST) see TODO
// // below
// template < typename TIVEL>
// extern Real _computeMaxIVel_OMP(Grid_t& grid);

template<typename TGrid>
Real FlowStep_LSRK3MPI<TGrid>::_computeMaxIVel()
{
    double maxIVel;

    // TODO: (fabianw@mavt.ethz.ch; Wed 06 Apr 2016 11:22:25 AM CEST) This is a
    // bit ugly due to different implementation strategy (LSRKstep struct) on
    // node layer.  I will change this once vectorized kernels are available
    // const double local_maxIVel = FlowStep_LSRK3::_computeMaxIVel();
    double local_maxIVel = _computeMaxIVel_OMP<MaxInterfaceVel_CPP_5eq>(grid);

    MPI_Comm mycart = grid.getCartComm();

    MPI_Allreduce(&local_maxIVel, &maxIVel, 1, MPI_DOUBLE, MPI_MAX, mycart);

    return maxIVel;
}

//Specialization
template<>
Real FlowStep_LSRK3MPI<GridMPI_t>::operator()(const Real max_dt, const Real current_time)
{
    set_constants();

    //here we just check stuff and compute the next dt
    if (verbosity && LSRK3data::step_id==0)
        cout << "Minimum grid spacing and smoothing length are: " << h_min << ", " << smoothlength << endl;

    LSRK3MPIdata::GSYNCH = parser("-gsync").asInt(omp_get_max_threads());

    MPI_Comm comm = grid.getCartComm();
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    const bool isroot = (0 == myrank);

    Timer timer;
    timer.start();
#ifdef _USE_HPM_
    if (LSRK3data::step_id>0) 	HPM_Start("dt");
#endif

    const Real maxSOS = _computeSOS();

    if (isnan(maxSOS))
    {
      if (isroot)
         cout << "Finishing... sos is NaN" << endl;

      abort();

      return 0;
    }
    if (maxSOS == 0)
    {
      if (isroot)
        cout << "Finishing... sos is 0" << endl;

      abort();

      return 0;
    }

    #ifdef _USE_HPM_
    if (LSRK3data::step_id>0) 		HPM_Stop("dt");
#endif
    const double t_sos = timer.stop();

    const Real dt_adv = CFL*h_min/maxSOS;
    Real dt = std::min(max_dt, dt_adv);

    const int mudt = parser("-mudt").asInt(1);
    Real dt_diff = 0.0;
    if (mudt && Simulation_Environment::MU_MAX>0) {
        dt_diff = (Real)(h_min*h_min/(12*Simulation_Environment::MU_MAX));
        dt = std::min(dt, dt_diff);
    }

    const int sigmadt = parser("-sigmadt").asInt(1);
    Real dt_surf = 0.0;
    if (sigmadt && Simulation_Environment::SIGMA>0.0)
    {
      // see Brackbill et al. JCP 1992
      const Real sumrho = Simulation_Environment::RHO1 + Simulation_Environment::RHO2;
      dt_surf = (Real)std::sqrt(sumrho*h_min*h_min*h_min/(4.0*M_PI*Simulation_Environment::SIGMA));
      dt = std::min(dt, dt_surf);
    }

    if (isroot)
    {
        cout << "sos max is " << maxSOS << "\n";
        cout << "advection dt is " << dt_adv << "\n";
        if (mudt && Simulation_Environment::MU_MAX>0)
          cout << "diffusion dt is " << dt_diff << "\n";
        if (sigmadt && Simulation_Environment::SIGMA>0.0)
          cout << "surface tension dt is " << dt_surf << "\n";
        if (parser("-source").asInt(0))
          cout << "source is active" << "\n";
        cout << "max dt is " << max_dt << "\n";
        cout << "dt is " << dt << "\n";
    }

    if (verbosity)
    {
        cout << "Dispatcher is " << LSRK3data::dispatcher << endl;
        cout << "Profiling information for sos is " << t_sos << endl;
    }

    if (isroot && maxSOS>1e6)
    {
        cout << "Speed of sound is too high. Is it realistic?" << endl;
        MPI_Abort(comm, 1);
    }

    if (isroot && dt <= std::numeric_limits<double>::epsilon())
    {
        cout << "Last time step encountered." << endl;

        return 0;
    }

    //now we perform an entire RK step
    const bool compact = parser("-sharp-form").asBool(true);

#ifdef _NONUNIFORM_BLOCK_

    if (parser("-kernels").asString("cpp")=="cpp")
        _process_LSRK3<Convection_CPP_HLLC_5eq_nonuniform, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_CPP, Update_State_CPP_5eq, InterfaceSharpening_CPP_cell_compact_5eq>(grid, dt, current_time);
// #if defined(_QPX_) || defined(_QPXEMU_)
//     else if (parser("-kernels").asString("cpp")=="qpx")
// #if defined(_BGQ_)
//         _process_LSRK3<Convection_QPX_HLLC_5eq, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_CPP, Update_State_CPP_5eq, InterfaceSharpening_CPP_cell_compact_5eq>(grid, dt, current_time);
// #else
//     _process_LSRK3<Convection_QPX_HLLC_5eq, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_QPX, Update_State_CPP_5eq, InterfaceSharpening_CPP_cell_compact_5eq>(grid, dt, current_time);
// #endif
// #endif
    else
    {
        cout << "combination not supported yet" << endl;
        MPI_Abort(comm, 1);
    }

#else

    if(compact)
    {
        if (parser("-kernels").asString("cpp")=="cpp")
            _process_LSRK3<Convection_CPP_HLLC_5eq, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_CPP, Update_State_CPP_5eq, InterfaceSharpening_CPP_cell_compact_5eq>(grid, dt, current_time);
#if defined(_QPX_) || defined(_QPXEMU_)
        else if (parser("-kernels").asString("cpp")=="qpx")
#if defined(_BGQ_)
            _process_LSRK3<Convection_QPX_HLLC_5eq, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_CPP, Update_State_CPP_5eq, InterfaceSharpening_CPP_cell_compact_5eq>(grid, dt, current_time);
#else
            _process_LSRK3<Convection_QPX_HLLC_5eq, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_QPX, Update_State_CPP_5eq, InterfaceSharpening_CPP_cell_compact_5eq>(grid, dt, current_time);
#endif
#endif
        else
        {
            cout << "combination not supported yet" << endl;
            MPI_Abort(comm, 1);
        }
    }
    else
    {
        if (parser("-kernels").asString("cpp")=="cpp")
            _process_LSRK3<Convection_CPP_HLLC_5eq, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_CPP, Update_State_CPP_5eq, InterfaceSharpening_CPP_5eq>(grid, dt, current_time);
#if defined(_QPX_) || defined(_QPXEMU_)
        else if (parser("-kernels").asString("cpp")=="qpx")
#if defined(_BGQ_)
            _process_LSRK3<Convection_QPX_HLLC_5eq, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_CPP, Update_State_CPP_5eq, InterfaceSharpening_CPP_5eq>(grid, dt, current_time);
#else
            _process_LSRK3<Convection_QPX_HLLC_5eq, Source_CPP, Diffusion_CPP, SurfaceTension_CPP, Update_QPX, Update_State_CPP_5eq, InterfaceSharpening_CPP_5eq>(grid, dt, current_time);
#endif
#endif
        else
        {
            cout << "combination not supported yet" << endl;
            MPI_Abort(comm, 1);
        }
    }

#endif /* _NONUNIFORM_BLOCK_ */


    LSRK3data::step_id++;

    return dt;
}
