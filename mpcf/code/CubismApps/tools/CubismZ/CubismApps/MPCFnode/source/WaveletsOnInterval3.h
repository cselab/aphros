/*
 *  WaveletsOnInterval3.h
 *  
 *
 *  Created by Diego Rossinelli on 3/1/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

namespace WaveletsOnInterval 
{	
#ifdef _FLOAT_PRECISION_
	typedef float FwtAp;
#else
	typedef double FwtAp;
#endif

template<bool lifting>
struct WI4
{		
    static inline FwtAp predict0_first(const FwtAp A0, const FwtAp A1, const FwtAp A2) 
	{
	    return 11./8 * A0 - 1./2 * A1 + 1./8 * A2;
	};

    static inline FwtAp predict1_first(const FwtAp A0, const FwtAp A1, const FwtAp A2)
	{
	    return 5./8 * A0 + 1./2 * A1 - 1./8 * A2; 
	};

    static inline FwtAp predict0_middle(const FwtAp A0, const FwtAp A1, const FwtAp A2) 
	{
	    return 1./8 * A0 + A1 - 1./8 * A2;
	};

    static inline FwtAp predict1_middle(const FwtAp A0, const FwtAp A1, const FwtAp A2)
	{
	    return -1./8 * A0 + A1 + 1./8 * A2;
	};

    static inline FwtAp predict0_last(const FwtAp A0, const FwtAp A1, const FwtAp A2) 
	{
	    return -1./8 * A0 + 1./2 * A1 + 5./8 * A2;
	};

    static inline FwtAp predict1_last(const FwtAp A0, const FwtAp A1, const FwtAp A2)
	{
	    return 1./8 * A0 - 1./2 * A1 + 11./8 * A2;
	};

    
    template<const int N, bool forward>
	static inline void transform(FwtAp data[N])
	{			
	    assert(N >= 8);
	    assert(N%2==0);
			
	    enum { Nhalf = N / 2 };
			
	    if (forward)
	    {
		FwtAp scalings[Nhalf];
		for(int i=0; i<Nhalf; i++)
		    scalings[i] = 0.5 * (data[2 * i] + data[2 * i + 1]);
		
		FwtAp details[Nhalf];
				
		// compute first detail
		details[0] = 0.5 * ((data[1] - data[0]) - (predict1_first(scalings[0], scalings[1], scalings[2]) - predict0_first(scalings[0], scalings[1], scalings[2])));
			
		// compute middle details
		for(int i = 1; i < Nhalf - 1; ++i)
		{
		    const int s = 2 * i;

		    details[i] = 0.5 * ((data[s + 1] - data[s]) - (predict1_middle(scalings[i - 1], scalings[i], scalings[i + 1]) - predict0_middle(scalings[i - 1], scalings[i], scalings[i + 1])));
		}
				
		// compute last detail
		details[Nhalf-1] = 0.5 * ((data[N-1] - data[N-2]) - (predict1_last(scalings[Nhalf-3], scalings[Nhalf-2], scalings[Nhalf-1]) - predict0_last(scalings[Nhalf-3], scalings[Nhalf-2], scalings[Nhalf-1])));
								
		copy(scalings, scalings + Nhalf, data);
		copy(details, details + Nhalf, data + Nhalf);
	    }
	    else
	    {
		FwtAp scalings[Nhalf], details[Nhalf];
		copy(data, data + Nhalf, scalings);
		copy(data + Nhalf, data + N, details);
								
		data[0] = predict0_first(scalings[0], scalings[1], scalings[2]) - details[0];
		data[1] = predict1_first(scalings[0], scalings[1], scalings[2]) + details[0];
						
		for(int i=1; i<Nhalf-1; i++)
		{
		    data[2 * i + 0] = predict0_middle(scalings[i-1], scalings[i], scalings[i+1]) - details[i];
		    data[2 * i + 1] = predict1_middle(scalings[i-1], scalings[i], scalings[i+1]) + details[i];
		}
		
		data[N-2] = predict0_last(scalings[Nhalf-3], scalings[Nhalf-2], scalings[Nhalf-1]) - details[Nhalf-1];
		data[N-1] = predict1_last(scalings[Nhalf-3], scalings[Nhalf-2], scalings[Nhalf-1]) + details[Nhalf-1];
	    }
	}
};

	template<typename WaveletType, int ROWSIZE, int COLSIZE>
	struct WaveletSweep
	{
		template<int BS, bool forward>
		inline void sweep1D(FwtAp data[BS][ROWSIZE])
		{
			for(int iy = 0; iy < BS; ++iy)
				WaveletType::template transform<BS, forward>(&data[iy][0]);
		}
		
		template<int BS>
		inline void xy_transpose(FwtAp data[BS][ROWSIZE])
		{
			for(int iy = 0; iy < BS; ++iy)
					for(int ix = iy + 1; ix < BS; ++ix)
					{
						const FwtAp temp = data[iy][ix];
						data[iy][ix] = data[ix][iy];
						data[ix][iy] = temp;
					}
		}
		
		template<int BS>
		inline void xz_transpose(FwtAp data[BS][COLSIZE][ROWSIZE])
		{
			for(int iy = 0; iy < BS; ++iy)	
				for(int iz = 0; iz < BS; ++iz)
					for(int ix = iz + 1; ix < BS; ++ix)
					{
						const FwtAp temp = data[iz][iy][ix];
						data[iz][iy][ix] = data[ix][iy][iz];
						data[ix][iy][iz] = temp;
					}
		}
		
		template<int BS, bool forward>
		inline void sweep2D(FwtAp data[BS][ROWSIZE])
		{			
			sweep1D<BS, forward>(data);
			xy_transpose<BS>(data);
			sweep1D<BS, forward>(data);
			//xy_transpose<BS>(data); we don't need this (cross the fingers)
		}
		
		template<int BS, bool bForward>
		inline void sweep3D(FwtAp data[BS][COLSIZE][ROWSIZE])
		{			
			if(bForward)
			{
				for(int iz = 0; iz < BS; ++iz)
					sweep2D<BS, true>(data[iz]);
				
				xz_transpose<BS>(data);
				
				for(int iz = 0; iz < BS; ++iz)
					for(int iy = 0; iy < BS; ++iy)
						WaveletType::template transform<BS, true>(&data[iz][iy][0]);
			}
			else
			{
				//xz_transpose(); we don't need this (cross the fingers)
				
				for(int iz = 0; iz < BS; ++iz)
					for(int iy = 0; iy < BS; ++iy)
						WaveletType::template transform<BS, false>(&data[iz][iy][0]);
				
				xz_transpose<BS>(data);
				
				for(int iz = 0; iz < BS; ++iz)
					sweep2D<BS, false>(data[iz]);
			}
		}
	};


}
