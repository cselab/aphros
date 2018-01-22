/*
 *  Test_SteadyState.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_STEADYSTATE_H_J5GZ8VHD
#define TEST_STEADYSTATE_H_J5GZ8VHD

#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>

#include "Tests.h"
#include "BlockInfo.h"
#include "Types.h"
#include "OutputProcessing.h"

#define __SERIALIZER_META_DATA(__SEP__) \
    restart_id __SEP__ t __SEP__ step_id __SEP__ dumper->m_tdump __SEP__ dumper->m_steplast __SEP__ dumper->m_timelast __SEP__ dumper->m_dumpcount


template <typename TGrid, typename TStepper, template <typename> class TSlice=Slice>
class Test_SteadyState: public Simulation
{
public:

    Test_SteadyState(ArgumentParser& P) :
        m_nonuniform(NULL), grid(NULL), stepper(NULL), dumper(NULL),
        isroot(true), _serializationKernel(NULL),
        parser(P), restart_id(0), t(0.0), step_id(0)
    { }

    virtual ~Test_SteadyState()
    {
        if (this->dumper)  delete dumper;
        if (this->stepper) delete stepper;
        if (this->grid)    delete grid;
        if (this->m_nonuniform)
                           delete m_nonuniform; // order matters here
    }

    virtual void run();
    virtual void setup();

    inline void verbosity(const int v) { VERBOSITY = v; }


protected:
    typedef void (*TSerializer)(const TGrid&, const int, const Real, const std::string, const std::string, const bool);
    typedef void (*TDeserializer)(TGrid&, const std::string, const std::string);

    // nonuniform grids
#ifdef _WENO3_
    NonUniformScheme<typename TGrid::BlockType, Weno3Coefficients>* m_nonuniform;
#else /* WENO5 */
    NonUniformScheme<typename TGrid::BlockType, Weno5Coefficients_Coralic>* m_nonuniform;
#endif /* _WENO3_ */

    // simulation data
    TGrid * grid;

    // time stepper
    TStepper * stepper;

    // output generation
    OutputProcessing<TGrid,TSlice> * dumper;

    // class state
    bool isroot;
    TSerializer _serializationKernel;

    // grid parameter
    int BPDX, BPDY, BPDZ;

    // stepper parameter
    int NSTEPS, SAVEPERIOD, ANALYSISPERIOD, VERBOSITY, REPORT_FREQ, REFRESHPERIOD;
    int LASTSAVE;
    Real CFL, TEND;

    // simulation parameter
    int restart_id;
    int MOLLFACTOR;
    bool bRESTART, bASCIIFILES, bAWK;
    bool bEXIT, bEXITSAVE;
    bool BC_PERIODIC[3];
    Real t, dt;
    int step_id;

    // helper
    ArgumentParser& parser;
    Profiler profiler;

    inline void _serialize_meta(std::ofstream& _of)
    {
        _of.setf(std::ios::scientific, std::ios::floatfield);
        _of << __SERIALIZER_META_DATA(<< '\t' <<);
    }
    inline void _serialize_meta(std::ifstream& _if)
    {
        _if >> __SERIALIZER_META_DATA(>>);
    }
    void _deserialize(TDeserializer kern);
    void _serialize(TSerializer kern);
    virtual void _post_save() {}
    virtual void _post_restart() {}
    virtual void _save()
    {
        if (LASTSAVE != step_id)
        {
            if (isroot)
            {
                std::cout << std::endl;
                this->_timestamp();
            }

            this->_serialize(_serializationKernel);
            LASTSAVE = step_id;

            this->_post_save();

            if (isroot)
                std::cout << std::endl;
        }
    }

    virtual void _setup_ic(TDeserializer _deserializer)
    {
        if (bRESTART)
        {
            this->_deserialize(_deserializer);
            this->_post_restart();
        }
        else
            this->_ic();
    }

    // user info
    void _infoBoard() const;

    inline void _timestamp()
    {
        time_t rawtime;
        std::time(&rawtime);
        struct tm* timeinfo = std::localtime(&rawtime);
        char buf[256];
        std::strftime(buf, 256, "%A, %h %d %Y, %r %Z %z", timeinfo);
        std::cout << buf << std::endl;
    }

    inline void _exitMsg(const std::string msg)
    {
        std::cout << "====> STOPPING SIMULATION..." << std::endl;
        std::cout << "====> Reason:....." << msg << std::endl;
        std::cout << "====> Step id:...." << step_id << " (of " << NSTEPS << ")" << std::endl;
        std::cout << "====> Time:......." << std::scientific << t << " (of " << std::scientific << TEND << ")" << std::endl;
        std::cout << "====> Next dt:...." << std::scientific << dt << std::endl;
        this->_timestamp();
        profiler.printSummary();
    }

    virtual void _setEnvironment()
    {
        parser.read_runtime_environment();

        REFRESHPERIOD  = parser.parseRuntime("refreshperiod").asInt();
        bEXIT          = parser.parseRuntime("exit").asBool();
        bEXITSAVE      = parser.parseRuntime("exitsave").asBool();

        NSTEPS         = parser.parseRuntime("nsteps").asInt();
        SAVEPERIOD     = parser.parseRuntime("saveperiod").asInt();
        ANALYSISPERIOD = parser.parseRuntime("analysisperiod").asInt();
        VERBOSITY      = parser.parseRuntime("verb").asInt();
        REPORT_FREQ    = parser.parseRuntime("report").asInt();

        CFL            = parser.parseRuntime("cfl").asDouble();
        TEND           = parser.parseRuntime("tend").asDouble();
        stepper->set_CFL(CFL);

        // OutputProcessing
        dumper->m_dumpperiod   = parser.parseRuntime("dumpperiod").asInt();
        dumper->m_dumpdt       = parser.parseRuntime("dumpdt").asDouble();
        dumper->m_bIO          = parser.parseRuntime("io").asBool();
        dumper->m_bVP          = parser.parseRuntime("vp").asBool();
        dumper->m_bHDF         = parser.parseRuntime("hdf").asBool();
        dumper->m_bHDF_SLICE   = parser.parseRuntime("hdf_slice").asBool();
        dumper->m_heavySkipStep= parser.parseRuntime("heavyskipstep").asInt();
        dumper->m_channels     = parser.parseRuntime("channels").asString();
    }

    // case specific
    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////         TEST STEADY STATE       ///////////////\n");
        printf("////////////////////////////////////////////////////////////\n");
        typedef typename TGrid::BlockType B;
        std::cout << "Domain size:   [" << BPDX*B::sizeX;
        std::cout << " x " << BPDY*B::sizeY;
        std::cout << " x " <<  BPDZ*B::sizeZ << "]" << std::endl;

        std::cout << "Domain extent: [" << Simulation_Environment::extents[0];
        std::cout << " x " << Simulation_Environment::extents[1];
        std::cout << " x " <<  Simulation_Environment::extents[2] << "]" << std::endl;
    }

    virtual void _setup_parameter();

    virtual void _init()
    {
        if (isroot)
        {
            // print run configuration to stdout
            _infoBoard();

            if (parser("dumplistall").asBool(false))
                dumper->list_registered();

            // write the initial set of parameter (everything seen up to this
            // point)
            parser.write_runtime_environment();
        }
    }
    virtual void _ic();

    virtual void _analysis()
    {
        // if (ANALYSISPERIOD != 0 && step_id%ANALYSISPERIOD == 0 && !bRESTART)
        // {
        //     profiler.push_start("ANALYSIS");
        //     profiler.pop_stop();
        // }
    }
    virtual void _pre_step() { }
    virtual void _post_step() { }
};

#undef __SERIALIZER_META_DATA


// static helpers
static bool exist_restart(const std::string suffix, const std::string& path)
{
    // check for existence of restart data for given suffix
    const std::string head = path + "/restart_HEAD";

    std::ifstream status(head.c_str());
    bool bstatus = status.good();
    if (bstatus)
    {
        int pair_id;
        status >> pair_id;
        status.close();
        std::ostringstream datafile;
        datafile << path << "/restart_data_" << pair_id << "." << suffix;
        std::ifstream datastatus(datafile.str().c_str());
        bstatus = datastatus.good();
        datastatus.close();
    }
    return bstatus;
}


// class implementation
template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::_deserialize(TDeserializer _kern)
{
    const std::string path = parser("-fpath").asString(".");
    const std::string meta = path + "/restart_HEAD";
    std::ifstream head(meta.c_str());
    if (!head.good())
    {
        if (isroot)
            std::cerr << "ERROR: Can not find restart_HEAD!" << std::endl;
        abort();
    }
    _serialize_meta(head);

    assert(restart_id >= 0);
    assert(t >= 0);
    assert(step_id >= 0);
    assert(dumper->m_tdump >= 0);
    assert(dumper->m_steplast >= 0);
    assert(dumper->m_timelast >= 0);

    const std::string suffix = parser("restart_format").asString("zbin");
    std::stringstream datafile;
    datafile << "restart_data_" << restart_id;
    if (_kern == NULL)
    {
        if (isroot)
            std::cerr << "ERROR: Can not find " << datafile.str() << "." << suffix << "!" << std::endl;
        abort();
    }
    else
        if (isroot)
        {
            std::cout << "DESERIALIZATION META DATA: restart_status_" << restart_id << std::endl;
            std::cout << "DESERIALIZATION: Time is " << std::scientific << t << " and step_id is " << step_id << std::endl;
        }
    _kern(*grid, datafile.str().c_str(), path.c_str());
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::_serialize(TSerializer _kern)
{
    const std::string path = parser("-fpath").asString(".");

    ++restart_id;

    std::stringstream datafile;
    datafile << "restart_data_" << restart_id;

    assert(_kern != NULL);
    _kern(*grid, step_id, t, datafile.str().c_str(), path.c_str(), false);

    if (isroot)
        std::cout << "SERIALIZATION: written file " << path << "/" << datafile.str() << std::endl;

    if (isroot)
    {
        std::stringstream statusfile;
        statusfile << path + "/restart_HEAD";
        std::ofstream status_0(statusfile.str().c_str());

        statusfile.str("");
        statusfile.clear();
        statusfile << path + "/restart_status_" << restart_id;
        std::ofstream status_1(statusfile.str().c_str());

        _serialize_meta(status_0); // HEAD
        _serialize_meta(status_1); // this one

        std::cout << "SERIALIZATION: written meta data " << statusfile.str() << std::endl;
    }

    // if specified, remove previous state
    if (isroot && restart_id > 1)
    {
        const bool removeold = parser("-removeold").asBool(false);
        const bool keep_meta = parser("-removeold_keep_meta").asBool(false);
        const std::string suffix = parser("restart_format").asString("zbin");
        if (removeold)
        {
            std::cout << "SERIALIZATION: removing previous state... " << std::endl;

            if (!keep_meta)
            {
                // remove old status file
                std::stringstream statusfile;
                statusfile << path + "/restart_status_" << restart_id - 1;
                std::remove(statusfile.str().c_str());
            }

            // remove old data file
            std::stringstream datafile;
            datafile << path + "/restart_data_" << restart_id - 1 << "." << suffix;
            std::remove(datafile.str().c_str());
        }
    }
    return;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::_infoBoard() const
{
    // print all run gonfiguration settings to stdout
    // (0) my header
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "-----------------------------  INFO BOARD  ---------------------------" << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    // (1) Runtime Arguments
    std::cout << "* COMPILETIME ARGUMENTS:" << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    {
        std::cout.width(50); std::cout.fill('.');
        std::cout << std::left << "Blocksize"  << ": " << _BLOCKSIZE_ << std::endl;
        std::cout.width(50); std::cout.fill('.');
#ifdef _FLOAT_PRECISION_
        std::cout << std::left << "Float Precision"  << ": SINGLE" << std::endl;
#else
        std::cout << std::left << "Float Precision"  << ": DOUBLE" << std::endl;
#endif /* _FLOAT_PRECISION_ */
        std::cout.width(50); std::cout.fill('.');
        std::cout << std::left << "Reciprocal Precision"  << ": " << _PREC_LEVEL_ << std::endl;
        std::cout.width(50); std::cout.fill('.');
        std::cout << std::left << "Code Fusion Level"  << ": " << _MICROFUSION_ << std::endl;
        std::cout.width(50); std::cout.fill('.');
#ifdef _BGQ_
        std::cout << std::left << "Architecture"  << ": PowerPC" << std::endl;
#else
        std::cout << std::left << "Architecture"  << ": x86" << std::endl;
#endif /* _BGQ_ */
        std::cout.width(50); std::cout.fill('.');
#if defined(_QPX_) || defined(_QPXEMU_)
        std::cout << std::left << "Vectorization Support"  << ": YES" << std::endl;
#else
        std::cout << std::left << "Vectorization Support"  << ": NO" << std::endl;
#endif /* defined(_QPX_) || defined(_QPXEMU_) */
        std::cout.width(50); std::cout.fill('.');
#ifdef _USE_HDF_
        std::cout << std::left << "HDF5 Support"  << ": YES" << std::endl;
#else
        std::cout << std::left << "HDF5 Support"  << ": NO" << std::endl;
#endif /* _USE_HDF_ */
        std::cout.width(50); std::cout.fill('.');
#ifdef _WENO3_
        std::cout << std::left << "WENO Scheme"  << ": WENO3" << std::endl;
#else
        std::cout << std::left << "WENO Scheme"  << ": WENO5" << std::endl;
#endif /* _WENO3_ */
        std::cout.width(50); std::cout.fill('.');
#ifdef _NOK_
        std::cout << std::left << "K*div(u) Enabled?"  << ": NO" << std::endl;
#else
        std::cout << std::left << "K*div(u) Enabled?"  << ": YES" << std::endl;
#endif /* _NOK_ */
    }
    // (2) The Lab
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "* THE LAB:" << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    {
        Lab dummy;
        std::cout.width(50); std::cout.fill('.');
        std::cout << std::left << "Lab Name"  << ": " << dummy.name() << std::endl;
        std::cout.width(50); std::cout.fill('.');
        std::cout << std::left << "x-Periodic"  << ": " << (dummy.is_xperiodic() ? "true" : "false" ) << std::endl;
        std::cout.width(50); std::cout.fill('.');
        std::cout << std::left << "y-Periodic"  << ": " << (dummy.is_yperiodic() ? "true" : "false" ) << std::endl;
        std::cout.width(50); std::cout.fill('.');
        std::cout << std::left << "z-Periodic"  << ": " << (dummy.is_zperiodic() ? "true" : "false" ) << std::endl;
        std::cout.width(50); std::cout.fill('.');
        std::cout << std::left << "Minimum grid spacing"  << ": " << std::scientific << stepper->get_hmin() << std::endl;
    }
    // (3) Runtime Arguments
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "* RUNTIME ARGUMENTS:" << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    parser.print_args();
    // (4) my footer
    std::cout << "----------------------------------------------------------------------" << std::endl;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::_setup_parameter()
{
    parser.mute();

    // required
    parser.set_strict_mode();
    BPDX       = parser("-bpdx").asInt();
    MOLLFACTOR = parser("-mollfactor").asInt();
    TEND       = parser("-tend").asDouble();
    CFL        = parser("-cfl").asDouble();
    parser.unset_strict_mode();

    // defaults
    BPDY           = parser("-bpdy").asInt(BPDX);
    BPDZ           = parser("-bpdz").asInt(BPDX);
    NSTEPS         = parser("-nsteps").asInt(0);
    bRESTART       = parser("-restart").asBool(false);

    VERBOSITY      = parser("-verb").asInt(0);
    REPORT_FREQ    = parser("-report").asInt(50);
    REFRESHPERIOD  = parser("-refreshperiod").asInt(15);
    bEXIT          = parser("-exit").asBool(false);
    bEXITSAVE      = parser("-exitsave").asBool(true);

    SAVEPERIOD     = parser("-saveperiod").asInt(0);
    LASTSAVE       = -1;
    ANALYSISPERIOD = parser("-analysisperiod").asInt(0);

    bAWK        = parser("-awk").asBool(false);
    bASCIIFILES = parser("-ascii").asBool(false);

    Simulation_Environment::RHO1   = parser("-rho1").asDouble(1.0);
    Simulation_Environment::RHO2   = parser("-rho2").asDouble(1.0);
    Simulation_Environment::P1     = parser("-p1").asDouble(1.0);
    Simulation_Environment::P2     = parser("-p2").asDouble(1.0);
    Simulation_Environment::GAMMA1 = parser("-g1").asDouble(1.4);
    Simulation_Environment::GAMMA2 = parser("-g2").asDouble(1.4);
    Simulation_Environment::PC1    = parser("-pc1").asDouble(0.0);
    Simulation_Environment::PC2    = parser("-pc2").asDouble(0.0);
    Simulation_Environment::C1     = std::sqrt(Simulation_Environment::GAMMA1*(Simulation_Environment::P1+Simulation_Environment::PC1)/Simulation_Environment::RHO1);
    Simulation_Environment::C2     = std::sqrt(Simulation_Environment::GAMMA2*(Simulation_Environment::P2+Simulation_Environment::PC2)/Simulation_Environment::RHO2);
    Simulation_Environment::extent = parser("-extent").asDouble(1.0);
    Simulation_Environment::MU1    = parser("-mu1").asDouble(0.0);
    Simulation_Environment::MU2    = parser("-mu2").asDouble(0.0);
    Simulation_Environment::SIGMA  = parser("-sigma").asDouble(0.0);

    if (parser("restart_format").asString("zbin") == "zbin") // we prioritize zbin
        _serializationKernel = &DumpZBin<TGrid, StreamerSerialization>;
#ifdef _USE_HDF_
    else if (parser("restart_format").asString("zbin") == "h5")
        _serializationKernel = &DumpHDF5<TGrid, StreamerSerialization>;
#endif
    else
    {
        if (isroot)
            std::cerr << "ERROR: No suitable restart format chosen!" << std::endl;
        abort();
    }
    assert(_serializationKernel != NULL);

    // some post computations
    {
        Lab dummy;

        BC_PERIODIC[0] = dummy.is_xperiodic();
        BC_PERIODIC[1] = dummy.is_yperiodic();
        BC_PERIODIC[2] = dummy.is_zperiodic();

        Simulation_Environment::BC_PERIODIC[0] = BC_PERIODIC[0];
        Simulation_Environment::BC_PERIODIC[1] = BC_PERIODIC[1];
        Simulation_Environment::BC_PERIODIC[2] = BC_PERIODIC[2];
    }

#ifdef _NONUNIFORM_BLOCK_
    Simulation_Environment::extents[0] = parser("-extentX").asDouble(1.0);
    Simulation_Environment::extents[1] = parser("-extentY").asDouble(1.0);
    Simulation_Environment::extents[2] = parser("-extentZ").asDouble(1.0);
#else
    const int BPD_MAX = max(max(BPDX, BPDY), BPDZ);
    Simulation_Environment::extents[0] = Simulation_Environment::extent * BPDX / (double) BPD_MAX;
    Simulation_Environment::extents[1] = Simulation_Environment::extent * BPDY / (double) BPD_MAX;
    Simulation_Environment::extents[2] = Simulation_Environment::extent * BPDZ / (double) BPD_MAX;
#endif /* _NONUNIFORM_BLOCK_ */

    Simulation_Environment::vol_bodyforce[0] = parser("-F1").asDouble(0.0);
    Simulation_Environment::vol_bodyforce[1] = parser("-F2").asDouble(0.0);
    Simulation_Environment::vol_bodyforce[2] = parser("-F3").asDouble(0.0);
    Simulation_Environment::grav_accel[0] = parser("-G1").asDouble(0.0);
    Simulation_Environment::grav_accel[1] = parser("-G2").asDouble(0.0);
    Simulation_Environment::grav_accel[2] = parser("-G3").asDouble(0.0);

    Simulation_Environment::MU_MAX = std::max(Simulation_Environment::MU1, Simulation_Environment::MU2);

    // some checks
    assert(TEND >= 0.0);
    assert(BPDX >= 1);
    assert(BPDY >= 1);
    assert(BPDZ >= 1);
    assert(CFL > 0 && CFL<1);
    assert(MOLLFACTOR > 0);
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::_ic()
{
    const double G1 = 1.0/(Simulation_Environment::GAMMA1-1);
    const double G2 = 1.0/(Simulation_Environment::GAMMA2-1);
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;

    std::vector<BlockInfo> vInfo = grid->getBlocksInfo();

    const double a1r1 = parser("a1r1").asDouble(0.5);
    const double a2r2 = parser("a2r2").asDouble(0.5);
    const double u    = parser("u").asDouble(0.0);
    const double v    = parser("v").asDouble(0.0);
    const double w    = parser("w").asDouble(0.0);
    const double p    = parser("p").asDouble(1.0);
    const double a2   = parser("a2").asDouble(0.0);
    const double a1   = 1.0 - a2;
    const double r    = a1r1 + a2r2;

    const double Gm = a1*G1 + a2*G2;
    const double Pm = a1*G1*F1 + a2*G2*F2;

    typedef typename TGrid::BlockType TBlock;

#pragma omp parallel for
    for(int i=0; i<(int)vInfo.size(); i++)
    {
        BlockInfo info = vInfo[i];
        TBlock& b = *(TBlock*)info.ptrBlock;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    b(ix, iy, iz).alpha1rho1= a1r1;
                    b(ix, iy, iz).alpha2rho2= a2r2;
                    b(ix, iy, iz).ru = r*u;
                    b(ix, iy, iz).rv = r*v;
                    b(ix, iy, iz).rw = r*w;
                    b(ix, iy, iz).energy = Gm*p + Pm + 0.5*r*(u*u + v*v + w*w);
                    b(ix, iy, iz).alpha2 = a2;
                }
    }
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::run()
{
    _init();

    dt = 1.0;
    while (true)
    {
        // concept:
        // (1) refresh runtime environment based on user input
        // (2) process output
        // (3) perform state analysis/statistics
        // (4) check for exit request
        // (5) check for serialization (there are _post_save() and
        // _post_restart() virtual methods for special needs in derived
        // classes)
        // (6) check for exit based on simulation state
        // (7) perform explicit time step (executes pre- and post-step methods
        // before and after if necessary for derived cases.  Override these
        // first, before overriding the run() method in derived classes)

        // (1)
        if (REFRESHPERIOD != 0 && step_id%REFRESHPERIOD == 0)
            _setEnvironment();

        // (2)
        const Real dtMax = (*dumper)(step_id, t, NSTEPS, TEND, profiler, true, bRESTART);

        // (3)
        _analysis();

        // (4)
        if (bEXIT)
        {
            if (SAVEPERIOD != 0 && bEXITSAVE) _save();
            if (isroot) _exitMsg("User request");
            break;
        }

        // (5)
        if (SAVEPERIOD != 0 && (step_id > 0 && step_id%SAVEPERIOD == 0) && !bRESTART)
        {
            profiler.push_start("SAVE");
            _save();
            profiler.pop_stop();
        }

        if (isroot && (REPORT_FREQ != 0 && step_id%REPORT_FREQ == 0))
        {
            _timestamp();
            profiler.printSummary();
        }

        // (6)
        if ((NSTEPS != 0 && step_id == NSTEPS) || dtMax <= std::numeric_limits<Real>::epsilon() || dt <= std::numeric_limits<Real>::epsilon())
        {
            if (SAVEPERIOD != 0 && bEXITSAVE) _save();
            if (isroot) _exitMsg("Simulation state");
            break;
        }

        if (isroot)
            std::cout << "--> Time is " << std::scientific << t << "; step_id is " << step_id << std::endl;

        // (7)
        _pre_step();

        profiler.push_start("STEP");
        dt = (*stepper)(dtMax, t);
        profiler.pop_stop();

        t += dt;
        ++step_id;
        bRESTART = false;

        _post_step();
    }
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::setup()
{
    _setup_parameter();

    if (isroot)
        _print_case_header();

#ifdef _NONUNIFORM_BLOCK_
    const double eX = Simulation_Environment::extents[0];
    const double eY = Simulation_Environment::extents[1];
    const double eZ = Simulation_Environment::extents[2];

#ifdef _WENO3_
    m_nonuniform = new NonUniformScheme<typename TGrid::BlockType, Weno3Coefficients>(0.0,eX,0.0,eY,0.0,eZ, BPDX, BPDY, BPDZ);
#else
    m_nonuniform = new NonUniformScheme<typename TGrid::BlockType, Weno5Coefficients_Coralic>(0.0,eX,0.0,eY,0.0,eZ, BPDX, BPDY, BPDZ);
#endif /* _WENO3_ */
    assert(m_nonuniform != NULL);

    // initialize scheme
    MeshDensityFactory mk(parser);
    m_nonuniform->init(mk.get_mesh_kernel(0), mk.get_mesh_kernel(1), mk.get_mesh_kernel(2));

    grid = new TGrid(&m_nonuniform->get_map_x(), &m_nonuniform->get_map_y(), &m_nonuniform->get_map_z(), BPDX, BPDY, BPDZ);
    assert(grid != NULL);

    vector<BlockInfo> infos = grid->getBlocksInfo();
    m_nonuniform->setup_coefficients(infos);

    m_nonuniform->print_mesh_statistics();

#else

    grid = new TGrid(BPDX, BPDY, BPDZ, Simulation_Environment::extent);
    assert(grid != NULL);
#endif /* _NONUNIFORM_BLOCK_ */

    stepper = new TStepper(*grid, CFL, parser, VERBOSITY);
    assert(stepper != NULL);

#ifdef _NONUNIFORM_BLOCK_
    stepper->set_hmin(m_nonuniform->minimum_cell_width());
#endif /* _NONUNIFORM_BLOCK_ */

    dumper = new OutputProcessing<TGrid,TSlice>(parser, *grid, isroot);
    assert(dumper != NULL);
    dumper->register_all(*grid);

    const std::string path = parser("-fpath").asString(".");
    TDeserializer _deserializer = NULL;
    if (exist_restart("zbin", path)) // we prioritize zbin (MPI topology must be the same as at serialization time)
        _deserializer = &ReadZBin<TGrid, StreamerSerialization>;
#ifdef _USE_HDF_
    else if (exist_restart("h5", path))
        _deserializer = &ReadHDF5<TGrid, StreamerSerialization>;
#endif

    _setup_ic(_deserializer);
}
#endif /* TEST_STEADYSTATE_H_J5GZ8VHD */
