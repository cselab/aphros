/*
 *  FullWaveletTransform.h
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#ifndef _FULLWAVELETTRANFORM_H_
#define _FULLWAVELETTRANFORM_H_ 1

#pragma once

#include <vector>
#include <algorithm>
#include <bitset>

#include "WaveletsOnInterval.h"	// 4th and 3rd order wavelets in C++ 
#if defined(_QPX_) || defined(_QPXEMU_)
#include "WaveletsOnIntervalQPX.h"	// 4th order wavelets (no lifting) in QPX  
#endif

using namespace std;

namespace WaveletsOnInterval 
{
	typedef WI4 ChosenWavelets;
	
	inline const char * _name (int wtype)
	{
		if (wtype == 0) return "None";
		if (wtype == 1) return "InterpWavelet4thOrder";
		if (wtype == 2) return "LiftedInterpWavelet4thOrder";
		if (wtype == 3) return "AverageInterpWavelet3rdOrder";
	}
	
	inline const char * ChosenWavelets_GetName(int wtype) { return _name(wtype); } 

	template<int BS, int ROWSIZE, int COLSIZE, int SLICESIZE>
#if 0	// peh: currently disabled,  #if defined(_QPX_) || defined(_QPXEMU_)
	struct FullTransformEngine : WaveletSweepQPX< ROWSIZE, COLSIZE>
#else
	struct FullTransformEngine : WaveletSweep< ChosenWavelets, ROWSIZE, COLSIZE>
#endif
	{		
		FullTransformEngine<BS/2, ROWSIZE, COLSIZE, SLICESIZE> child;
		
		inline void fwt(FwtAp data[SLICESIZE][COLSIZE][ROWSIZE], int wtype)
		{
			this->template sweep3D<BS, true>(data, wtype);
			
			child.fwt(data, wtype);
		}
		
		inline void iwt(FwtAp data[SLICESIZE][COLSIZE][ROWSIZE], int wtype)
		{
			child.iwt(data, wtype);
			
			this->template sweep3D<BS, false>(data, wtype);
		}

		template<typename DataType, int REFBS>
		int threshold(const FwtAp eps, bitset<REFBS * REFBS * REFBS>& mask_survivors, DataType * const buffer_survivors, const FwtAp data[SLICESIZE][COLSIZE][ROWSIZE])
		{
			enum { BSH = BS / 2 };
				
			const int survivors = child.template threshold<DataType, REFBS>(eps, mask_survivors, buffer_survivors, data);
			
			DataType * const buffer_start = buffer_survivors + survivors;
			int local_survivors = 0;

			for(int code = 1; code < 8; ++code)
			{
				const int xstart = BSH * (code & 1);
				const int ystart = BSH * (code / 2 & 1);
				const int zstart = BSH * (code / 4 & 1);
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
						{
							const int xsrc = xstart + ix;
							const int ysrc = ystart + iy;
							const int zsrc = zstart + iz;
														
							const FwtAp mydata = data[zsrc][ysrc][xsrc];
							
							const bool accepted = fabs(mydata) > eps; 
							
							const int dst = xsrc + REFBS * (ysrc + REFBS * zsrc);

							mask_survivors[dst] = accepted;
							
							if (accepted)
								buffer_start[local_survivors++] = (DataType)mydata;
						}
			}
	
			return survivors + local_survivors;			
		}
		
		
		template<typename DataType>
		void load(vector<DataType>& datastream, bitset<BS * BS * BS> mask, FwtAp data[SLICESIZE][COLSIZE][ROWSIZE])
		{			
			static const int BSH = BS / 2;
			
			for(int code = 7; code >= 1; --code)
			{
				const int xstart = BSH * (code & 1);
				const int ystart = BSH * (code / 2 & 1);
				const int zstart = BSH * (code / 4 & 1);
				
				for(int iz = BSH - 1; iz >= 0; --iz)
					for(int iy = BSH - 1; iy >= 0; --iy)
						for(int ix = BSH - 1; ix >= 0; --ix)
						{
							const int myx = xstart + ix;
							const int myy = ystart + iy;
							const int myz = zstart + iz;
							
							const int srcidx = myx + BS * (myy + BS * myz);
							
							assert(srcidx >= 0);
							assert(srcidx < BS * BS * BS);
							
							const bool eat = mask[srcidx];
							
							const DataType mydata = eat ? datastream.back() : 0;
							
							data[myz][myy][myx] = mydata;
							
							if (eat) 
								datastream.pop_back();						
						}
			}
			
			//code 0
			{
				bitset<BSH * BSH * BSH> childmask;
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
						{
							const int src = ix + BS * (iy + BS * iz);
							const int dst = ix + BSH * (iy + BSH * iz);
							
							childmask[dst] = mask[src];
						}
				
				child.load(datastream, childmask, data);
			}
		}
	};
	
	template<int ROWSIZE, int COLSIZE, int SLICESIZE>
	struct FullTransformEngine<4, ROWSIZE, COLSIZE, SLICESIZE>
	{
		enum { BS = 4 } ;
		
		void fwt(FwtAp data[SLICESIZE][COLSIZE][ROWSIZE], int wtype) { }
		
		void iwt(FwtAp data[SLICESIZE][COLSIZE][ROWSIZE], int wtype) { }
		
		template<typename DataType, int REFBS>
		int threshold(const FwtAp eps, bitset<REFBS * REFBS * REFBS>& mask_survivors, DataType * const buffer_survivors, const FwtAp data[SLICESIZE][COLSIZE][ROWSIZE])
		{	
			for(int iz = 0; iz < BS; ++iz)
				for(int iy = 0; iy < BS; ++iy)
					for(int ix = 0; ix < BS; ++ix)
						mask_survivors[ix + REFBS * (iy + REFBS * iz)] = true;
			
			for(int iz = 0, c = 0; iz < BS; ++iz)
				for(int iy = 0; iy < BS; ++iy)
					for(int ix = 0; ix < BS; ++ix, ++c)
						buffer_survivors[c] = data[iz][iy][ix];
			
			return BS * BS * BS;
		}
		
		template<typename DataType>
		void load(vector<DataType>& datastream, bitset<BS * BS * BS> mask, FwtAp data[SLICESIZE][COLSIZE][ROWSIZE])
		{
			assert(datastream.size() == BS * BS * BS);
			
			for(int iz = 0, c = 0; iz < BS; ++iz)
				for(int iy = 0; iy < BS; ++iy)
					for(int ix = 0; ix < BS; ++ix, ++c)
						data[iz][iy][ix] = datastream[c];
			
			datastream.clear();
		}
	};
	
	template<int BS>
	struct __attribute__((__aligned__(_ALIGNBYTES_))) FullTransform : FullTransformEngine<BS, BS, BS, BS>
	{
		FwtAp  data[BS][BS][BS];
		
		void fwt(int wtype)  { FullTransformEngine<BS, BS, BS, BS>::fwt(data, wtype); }
		
		void iwt(int wtype) { FullTransformEngine<BS, BS, BS, BS>::iwt(data, wtype); }
		
		template<typename DataType, int REFBS>
		int threshold(const FwtAp eps, bitset<REFBS * REFBS * REFBS>& mask_survivors, DataType * const buffer_survivors)
		{
			return FullTransformEngine<BS, BS, BS, BS>::template threshold<DataType, BS>(eps, mask_survivors, buffer_survivors, data);
		}

		template<typename DataType>
		void load(vector<DataType>& datastream, bitset<BS * BS * BS> mask)
		{
			FullTransformEngine<BS, BS, BS, BS>::template load<DataType>(datastream, mask, data);
		}
	};
}

#endif
