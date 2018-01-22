//
//  SurfaceTension_CPP.h
//  MPCFcore
//
//  Created by Babak Hejazialhosseini on 8/2/11.
//  Copyright 2011 ETH Zurich. All rights reserved.
//

#pragma once

#include "DivTensor_CPP.h"

class SurfaceTension_CPP: public virtual DivTensor_CPP
{
	struct AssumedType { Real a1r1, a2r2, ru, rv, rw, E, A2, dummy; };
	
public:
	SurfaceTension_CPP(const Real a, const Real dtinvh, 
			   const Real h, const Real sigma):
	DivTensor_CPP(a, dtinvh, h, sigma){ }
	
	//here we will "project" the phases
	void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
	
	void hpc_info(float& flop_convert, int& traffic_convert,
				  float& flop_corners, int& traffic_corner,
				  float& flop_tensorface, int& traffic_tensorface, 
				  float& flop_tdotu, int& traffic_tdotu,
				  float& flop_div, int& traffic_div,
				  float& flop_copyback, int& traffic_copyback,
				  size_t& footprint)
	{
		DivTensor_CPP::hpc_info(flop_convert, traffic_convert, flop_corners, traffic_corner, flop_tensorface, traffic_tensorface, 
					flop_tdotu, traffic_tdotu, flop_div, traffic_div, flop_copyback, traffic_copyback, footprint);
		
		const int ninputs = (int)powf(_BLOCKSIZE_ + 2, 3);
		
		flop_convert = 5 * ninputs;
		traffic_convert = (4 + 4) * sizeof(Real) * ninputs;
		
		footprint = sizeof(*this);
	}
};
