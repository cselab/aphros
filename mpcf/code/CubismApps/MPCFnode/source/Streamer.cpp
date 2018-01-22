/*
 *  Streamer.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 09/13/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#include "Streamer.h"

const std::string StreamerAlpha2::NAME = "Alpha2";
const std::string StreamerPressure::NAME = "Pressure";
const std::string StreamerDensity::NAME = "Density";
const std::string StreamerEnergy::NAME = "Energy";
const std::string StreamerQcriterion::NAME = "Qcriterion";
template <> const std::string StreamerVelocity<0>::NAME = "Velocity Ux";
template <> const std::string StreamerVelocity<1>::NAME = "Velocity Uy";
template <> const std::string StreamerVelocity<2>::NAME = "Velocity Uz";
const std::string StreamerVelocityVector::NAME = "Velocity vector";
const std::string StreamerVelocityMagnitude::NAME = "Velocity magnitude";
const std::string StreamerGradU::NAME = "Velocity gradient";
template <> const std::string StreamerVorticity<0>::NAME = "Vorticity Ox";
template <> const std::string StreamerVorticity<1>::NAME = "Vorticity Oy";
template <> const std::string StreamerVorticity<2>::NAME = "Vorticity Oz";
const std::string StreamerVorticityVector::NAME = "Vorticity vector";
const std::string StreamerVorticityMagnitude::NAME = "Vorticity magnitude";
const std::string StreamerAllConservative::NAME = "All conservative";
const std::string StreamerAllPrimitive::NAME = "All primitive";
const std::string StreamerK::NAME = "K";
const std::string StreamerDivU::NAME = "divU";
const std::string StreamerKDivU::NAME = "KdivU";
template <> const std::string StreamerAoScomponent<0>::NAME = "Component 0";
template <> const std::string StreamerAoScomponent<1>::NAME = "Component 1";
template <> const std::string StreamerAoScomponent<2>::NAME = "Component 2";
template <> const std::string StreamerAoScomponent<3>::NAME = "Component 3";
template <> const std::string StreamerAoScomponent<4>::NAME = "Component 4";
template <> const std::string StreamerAoScomponent<5>::NAME = "Component 5";
template <> const std::string StreamerAoScomponent<6>::NAME = "Component 6";
template <> const std::string StreamerAoScomponent<7>::NAME = "Component 7";
const std::string StreamerMach::NAME = "Mach";
const std::string StreamerSpeedOfSound::NAME = "Speed of sound";
const std::string StreamerSerialization::NAME = "Serialization";

const std::string StreamerAlpha2::EXT = "-a2";
const std::string StreamerPressure::EXT = "-p";
const std::string StreamerDensity::EXT = "-rho";
const std::string StreamerEnergy::EXT = "-E";
const std::string StreamerQcriterion::EXT = "-Qcrit";
template <> const std::string StreamerVelocity<0>::EXT = "-Ux";
template <> const std::string StreamerVelocity<1>::EXT = "-Uy";
template <> const std::string StreamerVelocity<2>::EXT = "-Uz";
const std::string StreamerVelocityVector::EXT = "-U";
const std::string StreamerVelocityMagnitude::EXT = "-IUI";
const std::string StreamerGradU::EXT = "-gradU";
template <> const std::string StreamerVorticity<0>::EXT = "-Omegax";
template <> const std::string StreamerVorticity<1>::EXT = "-Omegay";
template <> const std::string StreamerVorticity<2>::EXT = "-Omegaz";
const std::string StreamerVorticityVector::EXT = "-Omega";
const std::string StreamerVorticityMagnitude::EXT = "-IOmegaI";
const std::string StreamerAllConservative::EXT = "-Allc";
const std::string StreamerAllPrimitive::EXT = "-Allp";
const std::string StreamerK::EXT = "-K";
const std::string StreamerDivU::EXT = "-divU";
const std::string StreamerKDivU::EXT = "-KdivU";
template <> const std::string StreamerAoScomponent<0>::EXT = "-0";
template <> const std::string StreamerAoScomponent<1>::EXT = "-1";
template <> const std::string StreamerAoScomponent<2>::EXT = "-2";
template <> const std::string StreamerAoScomponent<3>::EXT = "-3";
template <> const std::string StreamerAoScomponent<4>::EXT = "-4";
template <> const std::string StreamerAoScomponent<5>::EXT = "-5";
template <> const std::string StreamerAoScomponent<6>::EXT = "-6";
template <> const std::string StreamerAoScomponent<7>::EXT = "-7";
const std::string StreamerMach::EXT = "-M";
const std::string StreamerSpeedOfSound::EXT = "-c";
const std::string StreamerSerialization::EXT = "";
