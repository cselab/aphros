/*
 *  Types.cpp
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/21/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include "Types.h"

Real Simulation_Environment::EPSILON;
Real Simulation_Environment::RHO1;
Real Simulation_Environment::RHO2;
Real Simulation_Environment::P1;
Real Simulation_Environment::P2;
Real Simulation_Environment::GAMMA1;
Real Simulation_Environment::GAMMA2;
Real Simulation_Environment::C1;
Real Simulation_Environment::C2;
Real Simulation_Environment::shock_pos;
Real Simulation_Environment::mach;
Real Simulation_Environment::reinit_tau_factor;
int Simulation_Environment::reinit_freq;
int Simulation_Environment::reinit_steps;
Real Simulation_Environment::PC1;
Real Simulation_Environment::PC2;
Real Simulation_Environment::extent;
Real Simulation_Environment::extents [3];
bool Simulation_Environment::BC_PERIODIC [3];
Real Simulation_Environment::vol_bodyforce[3] = {0};
Real Simulation_Environment::grav_accel[3] = {0};
Real Simulation_Environment::MU1;
Real Simulation_Environment::MU2;
Real Simulation_Environment::SIGMA;
Real Simulation_Environment::MU_MAX;

// Intel compiler does not allow in-class initialization of constant integal
// types.  Still, compilation is possible but linking fails? :(
const int FluidBlock::sizeX = _BLOCKSIZE_;
const int FluidBlock::sizeY = _BLOCKSIZE_;
const int FluidBlock::sizeZ = _BLOCKSIZE_;

const int FluidBlockNonUniform::sizeX = _BLOCKSIZE_;
const int FluidBlockNonUniform::sizeY = _BLOCKSIZE_;
const int FluidBlockNonUniform::sizeZ = _BLOCKSIZE_;


FluidElement operator * (const Real a, FluidElement gp)
{
    FluidElement out;

    out.alpha1rho1 = gp.alpha1rho1 * a;
    out.alpha2rho2 = gp.alpha2rho2 * a;
    out.ru         = gp.ru * a;
    out.rv         = gp.rv * a;
    out.rw         = gp.rw * a;
    out.energy     = gp.energy * a;
    out.alpha2     = gp.alpha2 * a;
    out.dummy      = gp.dummy * a;

    return out;
}

FluidElement operator + (FluidElement gpa, FluidElement gpb)
{
    FluidElement out;

    out.alpha1rho1 = gpa.alpha1rho1 + gpb.alpha1rho1;
    out.alpha2rho2 = gpa.alpha2rho2 + gpb.alpha2rho2;
    out.ru         = gpa.ru + gpb.ru;
    out.rv         = gpa.rv + gpb.rv;
    out.rw         = gpa.rw + gpb.rw;
    out.energy     = gpa.energy + gpb.energy;
    out.alpha2     = gpa.alpha2 + gpb.alpha2;
    out.dummy      = gpa.dummy  + gpb.dummy;

    return out;
}

FluidElement operator - (FluidElement gpa, FluidElement gpb)
{
    FluidElement out;

    out.alpha1rho1 = gpa.alpha1rho1 - gpb.alpha1rho1;
    out.alpha2rho2 = gpa.alpha2rho2 - gpb.alpha2rho2;
    out.ru         = gpa.ru - gpb.ru;
    out.rv         = gpa.rv - gpb.rv;
    out.rw         = gpa.rw - gpb.rw;
    out.energy     = gpa.energy - gpb.energy;
    out.alpha2     = gpa.alpha2 - gpb.alpha2;
    out.dummy      = gpa.dummy  - gpb.dummy;

    return out;
}
