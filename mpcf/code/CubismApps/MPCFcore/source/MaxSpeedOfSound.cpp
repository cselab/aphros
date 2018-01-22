/*
 *  MaxSpeedOfSound.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "common.h"
#include "MaxSpeedOfSound.h"

Real MaxSpeedOfSound_CPP::compute(const Real * const src, const int gptfloats) const
{
	const int N=_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	Real sos = 0;
	
	for(int i=0; i<N; i+=gptfloats)
	{
		const Real r = src[i];
		const Real u = src[i+1];
		const Real v = src[i+2];
		const Real w = src[i+3];
		const Real e = src[i+4];
		const Real G = src[i+5];
  		const Real P = src[i+6];

		assert(r>0);
		assert(e>0);

		assert(!isnan(r));
		assert(!isnan(u));
		assert(!isnan(v));
		assert(!isnan(w));
		assert(!isnan(e));
		assert(!isnan(G));
		assert(!isnan(P));

//		abort();
#if 1
		const Real p = (e - (u*u + v*v + w*w)*(0.5/r) - P)/G;

  		const Real c = sqrt(((p+P)/G+p)/r);
		
		assert(!isnan(p));
		assert(c > 0 && !isnan(c));

		sos = max(sos, c + max(max(abs(u), abs(v)), abs(w))/r);
#else
		const Real a = u*u + v*v + w*w;
		const Real ei = r*e - 0.5*a;

		const Real c = sqrt(ei*(1+G)-r*P)/(G*r);
		sos = max(sos, c + max(max(abs(u), abs(v)), abs(w))/r);

//		const Real c = sqrt(ei*(1+G)-r*P)/(G);
//		sos = max(sos, (c + max(max(abs(u), abs(v)), abs(w)))/r);
#endif
	}
    
	return sos;
}
