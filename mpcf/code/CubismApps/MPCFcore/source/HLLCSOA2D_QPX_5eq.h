/*
 *  HLLCSOA2D_QPX_5eq.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 03/03/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef HLLCSOA2D_QPX_5EQ_H_G9BFZP8H
#define HLLCSOA2D_QPX_5EQ_H_G9BFZP8H

#include "common.h"

class HLLCSOA2D_QPX_5eq
{
    const Real m_g1, m_g2;
    const Real m_pc1, m_pc2;

public:
    HLLCSOA2D_QPX_5eq(const Real g1, const Real g2, const Real pc1, const Real pc2) :
    m_g1(g1), m_g2(g2), m_pc1(pc1), m_pc2(pc2) {}

    void all(const TempSOA& r1minus, const TempSOA& r1plus,
             const TempSOA& r2minus, const TempSOA& r2plus,
             const TempSOA& vdminus, const TempSOA& vdplus,
             const TempSOA& v1minus, const TempSOA& v1plus,
             const TempSOA& v2minus, const TempSOA& v2plus,
             const TempSOA& pminus, const TempSOA& pplus,
             const TempSOA& A2minus, const TempSOA& A2plus,
             TempSOA& outam, TempSOA& outap, TempSOA& outas,
             TempSOA& outrho1, TempSOA& outrho2,
             TempSOA& outvd, TempSOA& outv1, TempSOA& outv2,
             TempSOA& oute,TempSOA& outA2) const;
};

#endif /* HLLCSOA2D_QPX_5EQ_H_G9BFZP8H */
