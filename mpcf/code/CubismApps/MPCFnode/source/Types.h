/*
 *  Types.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <fstream>
#include "math.h"
#include <iomanip>
#include <sstream>
#include <limits>
#include <vector>
#include <utility>
#include <iostream>
#include <omp.h>
#include <algorithm>
#ifdef _USE_NUMA_
#include <numa.h>
#endif

using namespace std;

#include <Grid.h>
#include <GridMorton.h>
#include <BlockLab.h>
#include <Profiler.h>
#include <ArgumentParser.h>
#include <SerializerIO.h>
#include <Timer.h>
#ifdef _USE_VTK_
#include <SerializerIO_ImageVTK.h>
#endif
#include <SerializerIO_VP.h>
#include <Timer.h>
#ifdef _USE_HDF_
#include <HDF5Dumper.h>
#endif
#include "NonUniform.h"
#include "GridTypes.h"

/* #define SETUP_MARKERS_IC \ */
/* const double mix_gamma = 1 + (G2*G1)/(G1*bubble+G2*(1-bubble)); \ */
/* const double mix_pinf  = (mix_gamma-1)/mix_gamma * (F1/G1*(1-bubble) + F2/G2*bubble); \ */
/* b(ix, iy, iz).G  = 1./(mix_gamma-1); \ */
/* b(ix, iy, iz).P = mix_gamma*mix_pinf/(mix_gamma-1); \ */
/* const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho; \ */
/* b(ix, iy, iz).energy   = pressure*b(ix, iy, iz).G + b(ix, iy, iz).P + ke; */

struct sort_pred {
    bool operator()(const std::pair<Real,Real> &left, const std::pair<Real,Real> &right) {
        return abs(left.first) < abs(right.first);
    }
};

class Simulation_Environment
{
public:
    static Real EPSILON, RHO1, RHO2, P1, P2, GAMMA1, GAMMA2, C1, C2;
    static Real shock_pos, mach;
    static Real reinit_tau_factor;
    static int reinit_freq, reinit_steps;
    static Real PC1, PC2;
    static Real extent;
    static Real extents [3];
    static bool BC_PERIODIC [3];
    static Real vol_bodyforce[3]; // volumetric body force (units Force/Volume)
    static Real grav_accel[3]; // gravitational acceleration (units Length/Time^2)
    static Real MU1, MU2, MU_MAX;
    static Real SIGMA;

    static Real pressure_5eq(const Real rho, const Real E, const Real ru, const Real rv, const Real rw, const Real alpha2)
    {
#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
        const Real a2 = max((Real)ALPHAEPS,min((Real)(1.0-ALPHAEPS),alpha2));
#else
        const Real a2 = alpha2;
#endif
        const Real a1 = 1.0 - a2;

        const Real g1m1Inv = 1.0/(GAMMA1-1.0);
        const Real g2m1Inv = 1.0/(GAMMA2-1.0);
        const Real gmix_m1 = 1.0/(a1*g1m1Inv + a2*g2m1Inv);
        const Real pcmix  = a1*GAMMA1*PC1*g1m1Inv + a2*GAMMA2*PC2*g2m1Inv;

#ifdef _CONVERTCLIP_
        const double actpres = gmix_m1 * (E - 0.5*(ru*ru + rv*rv + rw*rw)/rho) - gmix_m1*pcmix;
#ifdef _NOK_
        const Real pthresh = pcmix*gmix_m1/(gmix_m1+1.0);
#else
        // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 02:27:19 PM CEST) What
        // to do here if we use ALPHAEPS?
        const Real pthresh = a2 == static_cast<Real>(1.0) ? PC2 : (a2 == static_cast<Real>(0.0) ? PC1 : min(PC1,PC2));
#endif
        //const Real deltap = pthresh < (-static_cast<Real>(2.0)*actpres) ? (-static_cast<Real>(4.0)*actpres - pthresh) : static_cast<Real>(0.0);
        const Real deltap =  pthresh <= (-actpres + static_cast<Real>(PRESEPS)) ? (-actpres - pthresh + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);

        return (actpres + deltap);
#else
        return gmix_m1 * (E - 0.5*(ru*ru + rv*rv + rw*rw)/rho) - gmix_m1*pcmix;
#endif
    }

    static Real speedOfSound_5eq(const Real rho, const Real E, const Real ru, const Real rv, const Real rw, const Real a2)
    {
        const Real _g1 = Simulation_Environment::GAMMA1;
        const Real _g2 = Simulation_Environment::GAMMA2;
        const Real _pc1= Simulation_Environment::PC1;
        const Real _pc2= Simulation_Environment::PC2;
        const Real _g1m1Inv = 1.0/(_g1 - 1.0);
        const Real _g2m1Inv = 1.0/(_g2 - 1.0);

#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
        const Real alpha2 = std::max((Real)ALPHAEPS,std::min((Real)(1.0-ALPHAEPS),a2));
#else
        const Real alpha2 = a2;
#endif
        const Real alpha1 = static_cast<Real>(1.0) - alpha2;

        const Real rInv = static_cast<Real>(1.0) / rho;
        const Real GmixInv = static_cast<Real>(1.0) / (alpha1*_g1m1Inv + alpha2*_g2m1Inv);

#ifdef _CONVERTCLIP_
        const Real Pmix = std::max(static_cast<Real>(0.0),alpha1*_g1m1Inv*_g1*_pc1 + alpha2*_g2m1Inv*_g2*_pc2);
#else
        const Real Pmix = alpha1*_g1m1Inv*_g1*_pc1 + alpha2*_g2m1Inv*_g2*_pc2;
#endif /* _CONVERTCLIP_ */

        Real p = GmixInv*E - GmixInv*static_cast<Real>(0.5)*rInv*(ru*ru + rv*rv + rw*rw) - GmixInv*Pmix;

#ifdef _CONVERTCLIP_
#ifdef _NOK_
        const Real pthresh = Pmix*GmixInv/(GmixInv+1.0);
#else
        // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 02:46:04 PM CEST) What
        // to do here if ALPHAEPS is used?
        const Real pthresh = alpha2 == static_cast<Real>(1.0) ? _pc2 : (alpha2 == static_cast<Real>(0.0) ? _pc1 : std::min(_pc1, _pc2));
#endif
        const Real deltap = pthresh <= (-p + static_cast<Real>(PRESEPS)) ? (-p - pthresh + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
        p += deltap;
#endif /* _CONVERTCLIP_ */

#ifdef _NOK_
        return std::sqrt(((GmixInv + 1.0)*p + GmixInv*Pmix)*rInv);
#else

#if 1
        const Real denom1 = _g1*(p + _pc1);
        const Real denom2 = _g2*(p + _pc2);
        Real rc2 = (denom1 > 0) ? alpha1/denom1 : 0;
        rc2 += (denom2 > 0) ? alpha2/denom2 : 0;
        const Real c2 = rInv/rc2;
        return std::sqrt(c2);
#else
        return std::sqrt(rInv*_g1*_g2*(p+_pc1)*(p+_pc2) / (alpha1*_g2*(p+_pc2) + alpha2*_g1*(p+_pc1)));
#endif /* 0 */

#endif
    }

    static double heaviside(const double phi)
    {
        return (phi>0? 0:1);
    }

    static double heaviside_smooth(const double phi)
    {
        /*const double x = min((double)1, max((double)-1, phi*(((double)1)/EPSILON)));
        const double val_xneg = (((double)-0.5)*x - ((double)1))*x + ((double)0.5);
        const double val_xpos = (((double)+0.5)*x - ((double)1))*x + ((double)0.5);
        return (x<0 ? val_xneg : val_xpos);*/

        const double alpha = M_PI*min(1., max(0., 0.5*(phi+EPSILON)/EPSILON));
        return 0.5+0.5*cos(alpha);
    }

    static void getPostShockRatio(const Real pre_shock[3], const Real mach, const Real gamma, const Real pc, Real postShock[3])
    {
        const double Mpost = sqrt( (pow(mach,(Real)2.)*(gamma-1.)+2.) / (2.*gamma*pow(mach,(Real)2.)-(gamma-1.)) );
        postShock[0] = (gamma+1.)*pow(mach,(Real)2.)/( (gamma-1.)*pow(mach,(Real)2.)+2.)*pre_shock[0] ;
        postShock[2] = 1./(gamma+1.) * ( 2.*gamma*pow(mach,(Real)2.)-(gamma-1.))*pre_shock[2];
        const double preShockU = mach*sqrt(gamma*(pc+pre_shock[2])/pre_shock[0]);
        const double postShockU = Mpost*sqrt(gamma*(pc+postShock[2])/postShock[0]);
        postShock[1] = preShockU - postShockU;
    }

    /**
     * Computes post shock values given the post shock pressure
     *
     * @param array holding the pre shock values
     * @param pressure of the shock
     * @param array containing the computed post shock values
     *
     * @return void
     */
    static void getPostShockRatio(const Real pre_shock[3], const Real p_shock, Real postShock[3])
    {
        postShock[0] = pre_shock[0]*pow((double)p_shock/(3.31e4)+1,(1./7.));
        postShock[1] = pre_shock[1] + sqrt((postShock[0]-pre_shock[0])/(postShock[0]*pre_shock[0])*(p_shock-pre_shock[2]));
        postShock[2] = p_shock;
    }
};

class Simulation
{
public:
	virtual ~Simulation() {}
	virtual void run() = 0;
	virtual void setup() = 0;
};


struct StreamerGridPointASCII
{
	void operate(const FluidElement& input, ofstream& output) const
	{
        output << input.alpha1rho1 << " " << input.alpha2rho2 << " " << input.ru << " " << input.rv << " " <<
        input.rw << " " << input.energy << " " << input.alpha2 << endl;
	}

	void operate(ifstream& input, FluidElement& output) const
	{
		input >> output.alpha1rho1;
                input >> output.alpha2rho2;
		input >> output.ru;
		input >> output.rv;
		input >> output.rw;
		input >> output.energy;
		input >> output.alpha2;

		cout << "reading: " <<  output.alpha1rho1 << " "
                <<  output.alpha2rho2 << " "
		<<  output.ru << " "
		<<  output.rv << " "
		<<  output.rw << " "
		<<  output.energy << " "
                <<  output.alpha2 << "\n" ;
	}
};

struct StreamerGridPoint //dummy
{
    static const int channels = 7;

	void operate(const FluidElement& input, Real output[7]) const
	{
                const Real rho = input.alpha1rho1 + input.alpha2rho2;
                const Real alpha1 = 1.0 - input.alpha2;
                const Real gamma_m1 = 1.0/(alpha1/(Simulation_Environment::GAMMA1 - 1.0) + input.alpha2/(Simulation_Environment::GAMMA2 - 1.0));
                const Real gamma_pc_part = (alpha1 * Simulation_Environment::GAMMA1 * Simulation_Environment::PC1 / (Simulation_Environment::GAMMA1 - 1.0)
                                     + input.alpha2 * Simulation_Environment::GAMMA2 * Simulation_Environment::PC2 / (Simulation_Environment::GAMMA2 - 1.0));
		output[0] = rho;
		output[1] = input.ru/rho;
		output[2] = input.rv/rho;
		output[3] = input.rw/rho;
		output[4] = (input.energy-0.5*(input.ru*input.ru+input.rv*input.rv+input.rw*input.rw)/rho - gamma_pc_part) * gamma_m1;
		output[5] = input.energy; //TODO: URSULA just to have 7 values here (my by replaced by something else
  		output[6] = input.alpha2;
	}
};

template <> inline void FluidBlock::Write<StreamerGridPoint>(ofstream& output, StreamerGridPoint streamer) const
{
    output.write((const char *)&data[0][0][0], sizeof(FluidElement)*sizeX*sizeY*sizeZ);
}

template <> inline void FluidBlock::Read<StreamerGridPoint>(ifstream& input, StreamerGridPoint streamer)
{
    input.read((char *)&data[0][0][0], sizeof(FluidElement)*sizeX*sizeY*sizeZ);
}
