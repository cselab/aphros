/*
 *  Test_ChannelFlow.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger on 10/19/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_CHANNELFLOW_H
#define TEST_CHANNELFLOW_H

#include <vector>
#include "Test_SteadyState.h"

template <typename TGrid, typename TStepper, template <typename> class TSlice=Slice>
class Test_ChannelFlow: public virtual Test_SteadyState<TGrid,TStepper,TSlice>
{
public:
    Test_ChannelFlow(ArgumentParser& _parser) :
        Test_SteadyState<TGrid,TStepper,TSlice>(_parser) { }
    virtual ~Test_ChannelFlow() {}

protected:
    Real p0, rho0, u0, mu0, dpdx, delta, gamma, pc, mach;
    Real Re_tau, u_tau;

    typedef typename Test_SteadyState<TGrid,TStepper,TSlice>::TDeserializer TDeserializer;


    virtual void _setup_ic(TDeserializer _deserializer)
    {
        if (this->bRESTART)
        {
            this->_deserialize(_deserializer);
            this->_post_restart();

            // initialize forcing
            if (this->parser.check("-turbulent"))
                this->_set_forcing_turbulent();
            else
                this->_set_forcing_laminar();
        }
        else
            this->_ic();
    }

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////         TEST CHANNEL FLOW       ///////////////\n");
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
        this->_parameter_channel_flow();
    }

    virtual void _init()
    {
        Test_SteadyState<TGrid,TStepper,TSlice>::_init();
        this->_init_channel_flow();
    }

    virtual void _ic()
    {
        if (this->parser.check("-turbulent"))
        {
            this->_set_forcing_turbulent();
            this->_ic_channel_single_phase_turbulent();
        }
        else
        {
            this->_set_forcing_laminar();
            this->_ic_channel_single_phase();
        }
    }

    void _parameter_channel_flow()
    {
        // specific for this case
        gamma = this->parser("-gamma").asDouble(1.4);
        pc    = this->parser("-pc").asDouble(0.0);
        mach  = this->parser("-mach").asDouble(0.3);
        rho0  = this->parser("-rho0").asDouble(1.0);
        u0    = this->parser("-u0").asDouble(1.0);
        mu0   = this->parser("-mu0").asDouble(1.0);
        delta = this->parser("-delta").asDouble(1.0);
        if (this->parser.exist("-Re")) // laminar case
            Re_tau = this->parser("-Re").asDouble(50.0);
        else // turbulent case
            Re_tau = this->parser("-Re_tau").asDouble(180.0);

        // ensure single-phase setting
        Simulation_Environment::GAMMA1 = gamma;
        Simulation_Environment::GAMMA2 = gamma;
        Simulation_Environment::PC1 = pc;
        Simulation_Environment::PC2 = pc;

        Simulation_Environment::MU1 = mu0;
        Simulation_Environment::MU2 = 0.0;
        Simulation_Environment::MU_MAX = std::max(Simulation_Environment::MU1, Simulation_Environment::MU2);

        Simulation_Environment::vol_bodyforce[0] = 0.0;
        Simulation_Environment::vol_bodyforce[1] = 0.0;
        Simulation_Environment::vol_bodyforce[2] = 0.0;
    }

    void _init_channel_flow() const
    {
        if (this->isroot)
        {
            std::cout << "Mach   = " << this->mach << std::endl;
            std::cout << "p0     = " << this->p0 << std::endl;
            std::cout << "rho0   = " << this->rho0 << std::endl;
            std::cout << "u0     = " << this->u0 << std::endl;
            std::cout << "dpdx   = " << this->dpdx << std::endl;
            std::cout << "delta  = " << this->delta << std::endl;
            std::cout << "gamma  = " << Simulation_Environment::GAMMA1 << std::endl;
            std::cout << "pc     = " << Simulation_Environment::PC1 << std::endl;
            std::cout << "mu1    = " << Simulation_Environment::MU1 << std::endl;
            std::cout << "mu2    = " << Simulation_Environment::MU2 << std::endl;
            std::cout << "Re_tau = " << this->Re_tau << std::endl;

            std::cout << "Vol. Bodyforce X = " << Simulation_Environment::vol_bodyforce[0]  << std::endl;
            std::cout << "Vol. Bodyforce Y = " << Simulation_Environment::vol_bodyforce[1]  << std::endl;
            std::cout << "Vol. Bodyforce Z = " << Simulation_Environment::vol_bodyforce[2]  << std::endl;
        }
    }

    void _set_forcing_laminar();
    void _set_forcing_turbulent();
    void _match_velocity_energy(const double umeanIC);
    void _ic_channel_single_phase();
    void _ic_channel_single_phase_turbulent(const int offset=0);
};


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_ChannelFlow<TGrid,TStepper,TSlice>::_set_forcing_laminar()
{
    const double U_mean = 0.5*Re_tau*mu0/(rho0*delta); // Re_tau = Re_bulk here!
    dpdx = -3.0*mu0*U_mean/(delta*delta);
    Simulation_Environment::vol_bodyforce[0] = -dpdx; // set this pressure gradient as body force driver
    u0 = U_mean;

    // set pressure datum given Mach number and u0
    const double c_est = u0/mach;
    p0 = rho0/gamma*c_est*c_est;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_ChannelFlow<TGrid,TStepper,TSlice>::_set_forcing_turbulent()
{
    u_tau = Re_tau * mu0 / (rho0*delta);
    // \frac{\partial <p>}{\partial x} = \frac{d p_w}{d x} (c.f. Pope 7.1.2)
    dpdx = - rho0/delta * u_tau * u_tau; // (c.f. Pope 7.1.4)
    Simulation_Environment::vol_bodyforce[0] = -dpdx; // set this pressure gradient as body force driver

    // bulk Re
    const double Re = std::pow(Re_tau/0.09, 1.0/0.88); // (c.f. Pope 7.1.5)

    // mean bulk velocity
    u0 = 0.5 * Re * mu0 / (rho0*delta);

    // set pressure datum given Mach number and u0
    const double c_est = u0/mach;
    p0 = rho0/gamma*c_est*c_est;

    // if we restart from data with different Re_tau, adjust velocity and
    // energy (assumes same identical channel geometry, rho0 and mu0)
    if (this->bRESTART && this->parser.exist("-match_restart") && this->parser("-match_restart").asBool())
    {
        if (this->parser.exist("-match_re_tau"))
        {
            const double ic_Re_tau = this->parser("-match_re_tau").asDouble();
            const double ic_Re = std::pow(ic_Re_tau/0.09, 1.0/0.88); // (c.f. Pope 7.1.5)
            const double ic_u0 = 0.5 * ic_Re * mu0 / (rho0*delta);
            _match_velocity_energy(ic_u0);
        }
        else if (this->parser.exist("-match_mu0"))
        {
            const double ic_mu0 = this->parser("-match_mu0").asDouble();
            const double ic_u0 = 0.5 * Re * ic_mu0 / (rho0*delta);
            _match_velocity_energy(ic_u0);
        }
    }
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_ChannelFlow<TGrid,TStepper,TSlice>::_match_velocity_energy(const double umeanIC)
{
    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();

    const double G = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double F = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;

    typedef typename TGrid::BlockType TBlock;

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
                        // match the x-coordinate, leave fluctuations in y- and
                        // z-directions as is.
                        const double rho = b(ix, iy, iz).alpha1rho1 + b(ix, iy, iz).alpha2rho2;
                        b(ix, iy, iz).ru += (u0 - umeanIC)*rho;
                        b(ix, iy, iz).energy = (p0 + F)*G + 0.5 * (b(ix, iy, iz).ru * b(ix, iy, iz).ru + b(ix, iy, iz).rv * b(ix, iy, iz).rv + b(ix, iy, iz).rw * b(ix, iy, iz).rw)/rho;
                    }
        }
    }
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_ChannelFlow<TGrid,TStepper,TSlice>::_ic_channel_single_phase()
{
    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();

    const double G = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double F = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;

    typedef typename TGrid::BlockType TBlock;

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
                        Real p[3];
                        info.pos(p, ix, iy, iz);
                        b(ix, iy, iz).alpha1rho1 = rho0;
                        b(ix, iy, iz).alpha2rho2 = 0.0;

                        // parabolic (exact)
                        const double yd = (p[1]/delta - 1.0);
                        const double u = 1.5*u0*(1.0 - yd*yd);
                        b(ix, iy, iz).ru = u*rho0;

                        // uniform with equal mass flow rate as exact
                        // b(ix, iy, iz).ru = u0*rho0;

                        b(ix, iy, iz).rv = 0.0;
                        b(ix, iy, iz).rw = 0.0;
                        b(ix, iy, iz).energy = (p0 + F)*G + 0.5 * b(ix, iy, iz).ru * b(ix, iy, iz).ru/rho0;
                        b(ix, iy, iz).alpha2 = 0.0;
                    }
        }
    }
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_ChannelFlow<TGrid,TStepper,TSlice>::_ic_channel_single_phase_turbulent(const int offset)
{
    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();

    const double G = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double F = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;

    typedef typename TGrid::BlockType TBlock;

    // velocity fraction
    const double frac = this->parser("-U_fraction").asDouble(0.1);

    // init seed
    srand48(this->parser("-seed").asInt(42) + offset);

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
                        Real p[3];
                        info.pos(p, ix, iy, iz);
                        b(ix, iy, iz).alpha1rho1 = rho0;
                        b(ix, iy, iz).alpha2rho2 = 0.0;

                        const double yd = (p[1]/delta - 1.0);
                        const double u = 1.5*u0*(1.0 - yd*yd);
                        // const double u = u0;
                        // const double u = LEE_MOSER::LM_Umean_OneEightZero(p[1]/delta)*u_tau;

                        const double v = 0.0;
                        const double w = 0.0;

                        // add random fluctuations
                        // const double du = 0.0;
                        // const double dv = 0.0;
                        // const double dw = 0.0;
                        // const double du = frac*Umax_lam*drand48();
                        // const double du = frac*Umax_lam*(-1.0 + 2.0*drand48());
                        // const double dv = frac*Umax_lam*(-1.0 + 2.0*drand48());
                        // const double dw = frac*Umax_lam*(-1.0 + 2.0*drand48());

                        const double du = frac*u0*(-1.0 + 2.0*drand48());
                        const double dv = frac*u0*(-1.0 + 2.0*drand48());
                        const double dw = frac*u0*(-1.0 + 2.0*drand48());

                        const double U = u + du;
                        const double V = v + dv;
                        const double W = w + dw;
                        // const double V = v;
                        // const double W = w;

                        b(ix, iy, iz).ru = U*rho0;
                        b(ix, iy, iz).rv = V*rho0;
                        b(ix, iy, iz).rw = W*rho0;

                        b(ix, iy, iz).energy = (p0 + F)*G + 0.5 * rho0 * (U*U + V*V + W*W);
                        b(ix, iy, iz).alpha2 = 0.0;
                    }
        }
    }
}


namespace LEE_MOSER
{
    double LM_Umean_OneEightZero(double q)
    {
        const int Ndata = 96;
        const double yd[Ndata] = {
            1.110223024625157e-16,
            5.791461380333374e-05,
            2.141030191135096e-04,
            5.088966186284072e-04,
            9.825876412230539e-04,
            1.675417754204322e-03,
            2.627566689856309e-03,
            3.879140889489485e-03,
            5.412247554401284e-03,
            7.226453384460596e-03,
            9.321245632810338e-03,
            1.169603225078431e-02,
            1.435014205523744e-02,
            1.728282491824162e-02,
            2.049325197909280e-02,
            2.398051587857197e-02,
            2.774363101539001e-02,
            3.178153382474769e-02,
            3.609308307892944e-02,
            4.067706020984729e-02,
            4.553216965344320e-02,
            5.065703921585307e-02,
            5.605022046122787e-02,
            6.171018912110404e-02,
            6.763534552520467e-02,
            7.382401505355363e-02,
            8.027444860977029e-02,
            8.698482311541555e-02,
            9.395324202524546e-02,
            1.011777358632303e-01,
            1.086562627791852e-01,
            1.163867091258549e-01,
            1.243668900562930e-01,
            1.325945501413611e-01,
            1.410673640071792e-01,
            1.497829369923429e-01,
            1.587388058247237e-01,
            1.679324393176608e-01,
            1.773612390853474e-01,
            1.870225402772094e-01,
            1.969136123310693e-01,
            2.070316597448810e-01,
            2.173738228668202e-01,
            2.279371787035034e-01,
            2.387187417461112e-01,
            2.497154648141782e-01,
            2.609242399168145e-01,
            2.723418991311134e-01,
            2.839652154974980e-01,
            2.957909039317521e-01,
            3.078156221534808e-01,
            3.200359716307346e-01,
            3.324484985405327e-01,
            3.450496947450128e-01,
            3.578359987829319e-01,
            3.708037968762384e-01,
            3.839494239514293e-01,
            3.972691646754060e-01,
            4.107592545055342e-01,
            4.244158807536126e-01,
            4.382351836634480e-01,
            4.522132575017341e-01,
            4.663461516619241e-01,
            4.806298717807854e-01,
            4.950603808673232e-01,
            5.096336004437489e-01,
            5.243454116981767e-01,
            5.391916566487176e-01,
            5.541681393186459e-01,
            5.692706269223023e-01,
            5.844948510614016e-01,
            5.998365089314049e-01,
            6.152912645376163e-01,
            6.308547499206600e-01,
            6.465225663909910e-01,
            6.622902857720927e-01,
            6.781534516520057e-01,
            6.941075806428396e-01,
            7.101481636479068e-01,
            7.262706671361236e-01,
            7.424705344233169e-01,
            7.587431869600734e-01,
            7.750840256257712e-01,
            7.914884320284229e-01,
            8.079517698099662e-01,
            8.244693859566329e-01,
            8.410366121140243e-01,
            8.576487659065226e-01,
            8.743011522606657e-01,
            8.909890647321099e-01,
            9.077077868358066e-01,
            9.244525933790173e-01,
            9.412187517967892e-01,
            9.580015234895140e-01,
            9.747961651621934e-01,
            9.915979301650304e-01};
        const double Ud[Ndata] = {
            0.000000000000000e+00,
            1.054527346420741e-02,
            3.898150028900437e-02,
            9.264057361745806e-02,
            1.788297016728879e-01,
            3.048168764012739e-01,
            4.778105936207722e-01,
            7.049276244031449e-01,
            9.826535702095622e-01,
            1.310498745174586e+00,
            1.687704824580202e+00,
            2.113080599960335e+00,
            2.584774989317042e+00,
            3.100004373195065e+00,
            3.654779895855040e+00,
            4.243706300018299e+00,
            4.859931830114997e+00,
            5.495307682429419e+00,
            6.140767561712847e+00,
            6.786879590617578e+00,
            7.424476521405005e+00,
            8.045252641938248e+00,
            8.642230696990662e+00,
            9.210040607552013e+00,
            9.744998377122059e+00,
            1.024501391087453e+01,
            1.070938142498578e+01,
            1.113851203659422e+01,
            1.153365672387442e+01,
            1.189664875554680e+01,
            1.222968022110438e+01,
            1.253512087802303e+01,
            1.281538406627490e+01,
            1.307283702093434e+01,
            1.330974505927395e+01,
            1.352823955157148e+01,
            1.373030473089859e+01,
            1.391777480993175e+01,
            1.409232889493074e+01,
            1.425548107028632e+01,
            1.440857572773767e+01,
            1.455279522098163e+01,
            1.468917982624352e+01,
            1.481865160708460e+01,
            1.494202783963920e+01,
            1.506003393260543e+01,
            1.517331326548546e+01,
            1.528243123923937e+01,
            1.538788055948672e+01,
            1.549008706384810e+01,
            1.558942301410648e+01,
            1.568621884292019e+01,
            1.578077224609848e+01,
            1.587334286808345e+01,
            1.596412975734469e+01,
            1.605326477692437e+01,
            1.614083393295078e+01,
            1.622691086719516e+01,
            1.631157316543736e+01,
            1.639489659673339e+01,
            1.647695509123589e+01,
            1.655781587009242e+01,
            1.663752457029348e+01,
            1.671610261599261e+01,
            1.679355083037899e+01,
            1.686984536077985e+01,
            1.694494145794289e+01,
            1.701878494354969e+01,
            1.709131719124343e+01,
            1.716248304359947e+01,
            1.723223831568594e+01,
            1.730054835940958e+01,
            1.736737943380800e+01,
            1.743268741395911e+01,
            1.749640807217664e+01,
            1.755845313180477e+01,
            1.761871616470093e+01,
            1.767708464356820e+01,
            1.773344773526569e+01,
            1.778770118111745e+01,
            1.783974484741559e+01,
            1.788946716269731e+01,
            1.793674125922056e+01,
            1.798143758042809e+01,
            1.802342499001055e+01,
            1.806255782388764e+01,
            1.809867785606787e+01,
            1.813163300638841e+01,
            1.816129140022526e+01,
            1.818754492661940e+01,
            1.821029861421122e+01,
            1.822945651693691e+01,
            1.824491720935027e+01,
            1.825658665612607e+01,
            1.826439641260667e+01,
            1.826830831045079e+01};
        q = fabs(fabs(q-1.0) - 1.0);
        if (yd[Ndata-1] <= q) // extrapolate
            return Ud[Ndata-2] + (Ud[Ndata-1]-Ud[Ndata-2])/(yd[Ndata-1]-yd[Ndata-2])*(q-yd[Ndata-2]);
        else if (yd[0] > q) // extrapolate
            return Ud[0] + (Ud[1]-Ud[0])/(yd[1]-yd[0])*(q-yd[0]);
        else // interpolate
        {
            // q can not be zero due to cell centered grid representation
            const double qInv = 1.0/q;
            int k = 1;
            while (yd[k]*qInv < 1.0) ++k;
            return Ud[k-1] + (Ud[k]-Ud[k-1])/(yd[k]-yd[k-1])*(q-yd[k-1]);
        }
    }
}

#endif /*TEST_CHANNELFLOW_H*/
