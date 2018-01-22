/*
 *  Test_ChannelFlowMPI.h
 *  MPCFcluster
 *
 *  Created by Fabian Wermelinger on 10/20/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_CHANNELFLOWMPI_H_YVTWJDSI
#define TEST_CHANNELFLOWMPI_H_YVTWJDSI

#include "Test_ChannelFlow.h"
#include "Test_SteadyStateMPI.h"
#include "BlockProcessor_MPI.h"
#include "MPI_GridTransfer.h"
#include "StatisticsMPI_SPTurbulent.h"
#include "StatisticsMPI_SinglePhase.h"


template <typename TGrid, typename TStepper, template <typename> class TSlice=SliceMPI>
class Test_ChannelFlowMPI :
    public Test_ChannelFlow<TGrid,TStepper,TSlice>,
    public Test_SteadyStateMPI<TGrid,TStepper,TSlice>
{
public:
    Test_ChannelFlowMPI(const MPI_Comm comm, ArgumentParser& P) :
        Test_SteadyState<TGrid,TStepper,TSlice>(P),
        Test_ChannelFlow<TGrid,TStepper,TSlice>(P),
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>(comm,P)
    { }
    virtual ~Test_ChannelFlowMPI() {}


protected:
    bool m_bBodyforceByTime;
    int m_CHANNEL_BODYFORCE_UPDATE;
    int m_CHANNEL_STAT_UPDATE;
    int m_CHANNEL_STAT_DUMP;
    double m_Kp, m_Ki, m_Kd;

    // case specific
    virtual void _post_save()
    {
        if (this->isroot && m_bBodyforceByTime)
            _setBodyforce(true); // save bodyforce state too
    }

    virtual void _setEnvironment()
    {
        Test_SteadyState<TGrid,TStepper,TSlice>::_setEnvironment();

        m_bBodyforceByTime         = this->parser.parseRuntime("channel_temporal_bodyforce").asBool();
        m_CHANNEL_BODYFORCE_UPDATE = this->parser.parseRuntime("channel_bodyforce_update").asInt();
        m_CHANNEL_STAT_UPDATE      = this->parser.parseRuntime("channel_stat_update").asInt();
        m_CHANNEL_STAT_DUMP        = this->parser.parseRuntime("channel_stat_dump").asInt();

        m_Kp = this->parser.parseRuntime("Kp").asDouble();
        m_Ki = this->parser.parseRuntime("Ki").asDouble();
        m_Kd = this->parser.parseRuntime("Kd").asDouble();
    }

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////       TEST CHANNEL FLOW MPI     ///////////////\n");
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
        Test_ChannelFlow<TGrid,TStepper,TSlice>::_setup_parameter();
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_setup_parameter();

        m_bBodyforceByTime         = this->parser("channel_temporal_bodyforce").asBool(true);
        m_CHANNEL_BODYFORCE_UPDATE = this->parser("channel_bodyforce_update").asInt(1);
        m_CHANNEL_STAT_UPDATE      = this->parser("channel_stat_update").asInt(1);
        m_CHANNEL_STAT_DUMP        = this->parser("channel_stat_dump").asInt(1000);

        m_Kp = this->parser("Kp").asDouble(1.0);
        m_Ki = this->parser("Ki").asDouble(1.0);
        m_Kd = this->parser("Kd").asDouble(1.0);
    }

    virtual void _init()
    {
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_init();
        this->_init_channel_flow();

        m_hist = std::vector<SinglePhaseTurbulenceStatistics::AverageQoI>(this->BPDY*TGrid::BlockType::sizeY);
        m_channelTurnover = 0.0;

        // locations of planes for massflow computation
        const int Xall = this->XPESIZE * this->BPDX * TGrid::BlockType::sizeX;
        int locations[] = {0, Xall*0.25, Xall*0.5, Xall*0.75}; // evaluate at these locations
        m_loc = std::vector<int>(locations, locations + sizeof(locations)/sizeof(locations[0]));
    }

    virtual void _analysis()
    {
        if (this->ANALYSISPERIOD != 0 && this->step_id%this->ANALYSISPERIOD == 0 && !this->bRESTART)
        {
            this->profiler.push_start("ANALYSIS");
            SinglePhaseStatistics::dumpStatistics(*this->grid, this->step_id, this->t, this->dt, this->parser);
            this->profiler.pop_stop();
        }

        if (this->step_id % m_CHANNEL_STAT_UPDATE == 0 && this->step_id != 0)
        {
            this->profiler.push_start("CHANNEL STATISTICS UPDATE");
            SinglePhaseTurbulenceStatistics::updateStatistics(*this->grid, m_hist);
            this->profiler.pop_stop();
        }

        if (this->step_id % m_CHANNEL_STAT_DUMP == 0 && this->step_id != 0)
        {
            this->profiler.push_start("CHANNEL STATISTICS DUMP");
            SinglePhaseTurbulenceStatistics::dumpStatistics(*this->grid, m_hist,
                    this->step_id, this->t, this->dt, this->parser);
            this->profiler.pop_stop();
        }
    }

    virtual void _pre_step()
    {
        // compute current massflow per unit z-direction
        this->profiler.push_start("COMPUTE MASSFLOW RATE");
        m_mdot = this->_computeFlowRate();
        this->profiler.pop_stop();

        if (m_bBodyforceByTime && this->step_id % m_CHANNEL_BODYFORCE_UPDATE == 0)
        {
            this->profiler.push_start("CHANNEL BODYFORCE UPDATE");
            _setBodyforce();
            this->profiler.pop_stop();
        }
    }

    virtual void _post_step()
    {
        // TODO: (fabianw@mavt.ethz.ch; Mon 08 Aug 2016 02:47:53 PM CEST) This
        // assumes incompressible flow
        m_channelTurnover += m_mdot * this->dt/(2.0 * this->rho0 * this->delta * Simulation_Environment::extents[0]);
        if (this->isroot) printf("channel turnover = %f\n", m_channelTurnover);
    }

    double _computeFlowRate();
    void _setBodyforce(const bool bSave=false);
    void _interpolateGrid();

private:
    Real m_mdot;
    std::vector<SinglePhaseTurbulenceStatistics::AverageQoI> m_hist;
    std::vector<int> m_loc;
    double m_channelTurnover;
};


template <typename TGrid, typename TStepper, template <typename> class TSlice>
double Test_ChannelFlowMPI<TGrid,TStepper,TSlice>::_computeFlowRate()
{
    typedef typename TGrid::BlockType TBlock;
    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();

    const double h  = vInfo[0].h_gridpoint;
    const double LzInv = 1.0/(this->ZPESIZE * this->BPDZ * TBlock::sizeZ);

    int peidx[3];
    this->grid->peindex(peidx);

    std::vector<int> myloc;
    const int gstart = peidx[0] * this->BPDX * TBlock::sizeX;
    const int gend = gstart + this->BPDX * TBlock::sizeX;
    for (size_t i = 0; i < m_loc.size(); ++i)
        if (m_loc[i] >= gstart && m_loc[i] < gend)
            myloc.push_back(m_loc[i]);

    std::vector<double> Usum(myloc.size(), 0.0);

    double t0, t1;

    t0 = omp_get_wtime();
    for (size_t l = 0; l < myloc.size(); l++)
    {
        const int gx = myloc[l];
        double& sum = Usum[l];

        // TODO: (fabianw@mavt.ethz.ch; Fri 21 Oct 2016 02:53:59 PM CEST) use
        // OMP here
        for(int i=0; i<vInfo.size(); i++)
        {
            const BlockInfo& info = vInfo[i];
            TBlock& b = *(TBlock *)info.ptrBlock;

            const int bstart = vInfo[i].index[0]*TBlock::sizeX;
            const int bend = bstart + TBlock::sizeX;
            if (gx < bstart || bend <= gx) continue;

            const int ix = (gx % (this->BPDX*TBlock::sizeX)) % TBlock::sizeX;

            for(int iz=0; iz<TBlock::sizeZ; iz++)
                for(int iy=0; iy<TBlock::sizeY; iy++)
                    sum += b(ix,iy,iz).ru;
        }
    }
    t1 = omp_get_wtime();
    if (this->isroot) printf("massflow: A took %lf seconds\n", t1-t0);

    MPI_Barrier(this->m_comm_world);

    t0 = omp_get_wtime();

    int worldRank, worldSize;
    MPI_Comm_size(this->m_comm_world, &worldSize);
    MPI_Comm_rank(this->m_comm_world, &worldRank);

    MPI_Comm comm_subgroup;
    MPI_Comm_split(this->m_comm_world, peidx[0], worldRank, &comm_subgroup);

    int subrank;
    MPI_Comm_rank(comm_subgroup, &subrank);

    int myCount = 0;
    if (!Usum.empty())
    {
        // reduction on comm_subgroups
        if (0==subrank)
        {
            myCount = Usum.size();
            MPI_Reduce(MPI_IN_PLACE, Usum.data(), Usum.size(), MPI_DOUBLE, MPI_SUM, 0, comm_subgroup);
        }
        else
            MPI_Reduce(Usum.data(), Usum.data(), Usum.size(), MPI_DOUBLE, MPI_SUM, 0, comm_subgroup);
    }

    std::vector<int> howMany(worldSize, 0);
    std::vector<int> offsets(worldSize, 0);
    MPI_Gather(&myCount, 1, MPI_INT, howMany.data(), 1, MPI_INT, 0, this->m_comm_world);
    for (size_t i = 1; i < offsets.size(); ++i)
        offsets[i] = offsets[i-1] + howMany[i-1];

    std::vector<double> gUsum(m_loc.size(),0.0);
    MPI_Gatherv(Usum.data(), myCount, MPI_DOUBLE, gUsum.data(), howMany.data(), offsets.data(), MPI_DOUBLE, 0, this->m_comm_world);

    MPI_Comm_free(&comm_subgroup);

    double mdot_avgOverLoc = 0.0;
    if (this->isroot)
    {
        std::ofstream out("massflow.dat", std::ios::app);
        out << this->step_id;
        out.setf(std::ios::scientific, std::ios::floatfield);
        out.precision(12);
        out << '\t' << this->t;
        for (size_t i = 0; i < gUsum.size(); ++i)
        {
            mdot_avgOverLoc += gUsum[i]*LzInv*h;
            out << '\t' << gUsum[i]*LzInv*h;
        }
        out << std::endl;
        out.close();
        mdot_avgOverLoc /= gUsum.size();
    }

    MPI_Bcast(&mdot_avgOverLoc, 1, MPI_DOUBLE, 0, this->m_comm_world);

    t1 = omp_get_wtime();
    if (this->isroot) printf("massflow: B took %lf seconds\n", t1-t0);

    return mdot_avgOverLoc;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_ChannelFlowMPI<TGrid,TStepper,TSlice>::_setBodyforce(const bool bSave)
{
    static double sum = 0.0;
    static double old = 0.0;
    static size_t invocation = 0;
    if (bSave)
    {
        ofstream savefile("bodyforce.par");
        savefile << setprecision(16) << sum << std::endl;
        savefile << setprecision(16) << old << std::endl;
        savefile << setprecision(16) << this->dt << std::endl;
        savefile << setprecision(16) << Simulation_Environment::vol_bodyforce[0] << std::endl;
        savefile << invocation << std::endl;
        savefile << "Iteration = " << this->step_id << std::endl;
        savefile.close();
        return;
    }
    if (invocation == 0)
    {
        ifstream loadfile("bodyforce.par");
        if (loadfile.good())
        {
            loadfile >> sum;
            loadfile >> old;
            loadfile >> this->dt;
            loadfile >> Simulation_Environment::vol_bodyforce[0];
            loadfile >> invocation;
        }
    }

    // compute error (deviation from target value)
    const double delta = this->delta;
    const double e = 1.5 * this->mu0 / (this->rho0 * delta*delta*delta)*(2.0 * this->rho0 * delta * this->u0 - m_mdot);

    // PID contributions (parameter are set online)
    if (invocation > 0)
    {
        // integral
        sum += 0.5*m_Ki*(old + e) * this->dt;
        // differential
        const double dterm = m_Kd*(e - old) / this->dt;
        // controller output
        const double u = m_Kp*e + sum + dterm;

        // update control variable
        Simulation_Environment::vol_bodyforce[0] += u;
    }
    else
        old = e;
    ++invocation;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_ChannelFlowMPI<TGrid,TStepper,TSlice>::_interpolateGrid()
{
    // construct (previous) coarser grid
    // we assume that the resolution increases by a factor of 2
    assert(this->BPDX%2 == 0 && this->BPDY%2 == 0 && this->BPDZ%2 == 0);
    const int coarse_bpdx = this->BPDX/2;
    const int coarse_bpdy = this->BPDY/2;
    const int coarse_bpdz = this->BPDZ/2;

    // you must read the coarse grid using the MPI topology of the fine grid!
    // Thus, you can not use ZBIN.
    TGrid * const coarse_grid = new TGrid(this->XPESIZE, this->YPESIZE, this->ZPESIZE,
            coarse_bpdx, coarse_bpdy, coarse_bpdz,
            Simulation_Environment::extent, this->m_comm_world);

#ifdef _USE_HDF_
    const std::string path = this->parser("-fpath").asString(".");
    this->parser.set_strict_mode();
    const std::string fNameCoarse = this->parser("-coarse_data").asString();
    this->parser.unset_strict_mode();

    ReadHDF5_MPI<TGrid, StreamerSerialization>(*coarse_grid, fNameCoarse, path);
    DumpHDF5_MPI<TGrid, StreamerAllPrimitive>(*coarse_grid, 0, 0.0, "coarseRead");
#else
    if (this->isroot)
        std::cout << "Enable HDF for grid transfer. Abort now ..." << std::endl;
    abort();
#endif

    grid_transfer<TGrid> coarse_to_fine(*this->grid,
            coarse_bpdx, coarse_bpdy, coarse_bpdz,
            this->BPDX, this->BPDY, this->BPDZ);

    // This requires the that coarse block on THIS process contains the fine
    // blocks
    //
    // o------x------o
    // |      |      |
    // |      |      |
    // x------x------x
    // |      |      |
    // |      |      |
    // o------x------o
    // Coarse block corners "o", fine blocks corners "x", both must be in the
    // memory space of this process
    process< LabMPI >(coarse_to_fine, *coarse_grid, 0, 0);

    // clean up
    delete coarse_grid;

    return;
}

#endif /* TEST_CHANNELFLOWMPI_H_YVTWJDSI */
