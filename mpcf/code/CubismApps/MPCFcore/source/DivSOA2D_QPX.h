/*
 *  DivSOA2D_QPX.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 03/03/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef DIVSOA2D_QPX_H_RGAODOCK
#define DIVSOA2D_QPX_H_RGAODOCK

#include "common.h"

class DivSOA2D_QPX
{
public:
    void xrhs(const TempSOA& flux, OutputSOA& rhs) const;
    void yrhs(const TempSOA& flux, OutputSOA& rhs) const;
    void zrhs(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs) const;
};


class DivSOA2D_HLLC_QPX_5eq
{
    const Real m_g1, m_g2;
    const Real m_pc1, m_pc2;

public:
    DivSOA2D_HLLC_QPX_5eq(const Real g1, const Real g2, const Real pc1, const Real pc2) :
    m_g1(g1), m_g2(g2), m_pc1(pc1), m_pc2(pc2) {}

    void xextraterm(const TempSOA& um, const TempSOA& up,
                    const TempSOA& A2m, const TempSOA& A2p,
                    const TempSOA& pm, const TempSOA& pp,
                    const TempSOA& am, const TempSOA& ap, const TempSOA& as,
                    OutputSOA& divu, OutputSOA& sumA2, OutputSOA& sumK);

    void yextraterm(const TempSOA& vm, const TempSOA& vp,
                    const TempSOA& A2m, const TempSOA& A2p,
                    const TempSOA& pm, const TempSOA& pp,
                    const TempSOA& am, const TempSOA& ap, const TempSOA& as,
                    OutputSOA& divu, OutputSOA& sumA2, OutputSOA& sumK);

    void zextraterm(const TempSOA& wm0, const TempSOA& wp0, const TempSOA& wm1, const TempSOA& wp1,
                    const TempSOA& A2m, const TempSOA& A2p, const TempSOA& pm, const TempSOA& pp,
                    const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1,
                    const TempSOA& as0, const TempSOA& as1,
                    OutputSOA& divu, OutputSOA& sumA2, OutputSOA& sumK);
};

#endif /* DIVSOA2D_QPX_H_RGAODOCK */
