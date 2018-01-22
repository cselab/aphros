/*
 *  Update_State_CPP_5eq.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 03/19/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef UPDATE_STATE_CPP_5EQ_H_EWXH62ZF
#define UPDATE_STATE_CPP_5EQ_H_EWXH62ZF
#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "common.h"

class Update_State_CPP_5eq
{
protected:

    const Real g1, g2, pc1, pc2;
    const Real g1m1Inv, g2m1Inv;
    const Real m_alpha;
    const Real m_beta;

public:

    Update_State_CPP_5eq(const Real gamma1, const Real gamma2, const Real pc1, const Real pc2, const Real alpha=2.0, const Real beta=4.0): g1(gamma1), g2(gamma2), pc1(pc1), pc2(pc2), g1m1Inv(1.0/(g1-1.0)), g2m1Inv(1.0/(g2-1.0)), m_alpha(alpha), m_beta(beta) {}

    void compute(Real * const dst, const int gptfloats) const;
};

#endif /* UPDATE_STATE_CPP_5EQ_H_EWXH62ZF */
