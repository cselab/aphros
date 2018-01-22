/*
 *  MaxSpeedOfSound_CPP_5eq.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 10/22/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef MAXSPEEDOFSOUND_CPP_5EQ_H_IWC2DQZ9
#define MAXSPEEDOFSOUND_CPP_5EQ_H_IWC2DQZ9

#include "MaxSpeedOfSound.h"

class MaxSpeedOfSound_CPP_5eq : public MaxSpeedOfSound_CPP
{
public:
    MaxSpeedOfSound_CPP_5eq(const Real g1=1.4, const Real g2=1.4, const Real pc1=0.0, const Real pc2=0.0) :
        MaxSpeedOfSound_CPP(), _g1(g1), _g2(g2), _pc1(pc1), _pc2(pc2), _g1m1Inv(1.0/(g1-1.0)), _g2m1Inv(1.0/(g2-1.0)) {}

    virtual Real compute(const Real * const src, const int gptfloats) const;

protected:
    const Real _g1, _g2;
    const Real _pc1, _pc2;
    const Real _g1m1Inv, _g2m1Inv; // inverse of constant (gamma_k - 1) for k = {1,2}
};

#endif /* MAXSPEEDOFSOUND_CPP_5EQ_H_IWC2DQZ9 */
