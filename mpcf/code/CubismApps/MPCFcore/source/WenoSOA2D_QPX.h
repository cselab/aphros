/*
 **  WenoSOA2D_QPX.h
 **
 **
 **  Created by Diego Rossinelli on 2/12/13.
 **  Copyright 2013 ETH Zurich. All rights reserved.
 **
 **/
#ifndef WENOSOA2D_QPX_H_67YJWBOH
#define WENOSOA2D_QPX_H_67YJWBOH

#include "common.h"

#include "../../MPCFthread/source/Weno_CPP.h"
#include "../../MPCFthread/source/Weno_QPX.h"

#ifndef _WENO3_
#ifdef _ACCURATEWENO_
typedef WenoQPX_MinusFunctor SoupMinus;
typedef WenoQPX_PlusFunctor SoupPlus;
#else
typedef WenoQPX_MinusPanos SoupMinus;
typedef WenoQPX_PlusPanos SoupPlus;
#endif
#else
#include "../../MPCFthread/source/Weno_QPX_3rdOrder.h"
typedef WenoQPX_3rdOrder_Minus SoupMinus;
typedef WenoQPX_3rdOrder_Plus SoupPlus;
#endif

class WenoSOA2D_QPX
{
public:

	void xcompute(const InputSOA& in, TempSOA& outm, TempSOA& outp) const;
	void ycompute(const InputSOA& in, TempSOA& outm, TempSOA& outp) const;
	void zcompute(const int r, const RingInputSOA& in, TempSOA& outm, TempSOA& outp) const;

    void xcompute_P(const InputSOA& in, TempSOA& outm, TempSOA& outp) const;
	void ycompute_P(const InputSOA& in, TempSOA& outm, TempSOA& outp) const;
	void zcompute_P(const int r, const RingInputSOA& in, TempSOA& outm, TempSOA& outp) const;
};

#endif /* WENOSOA2D_QPX_H_67YJWBOH */
