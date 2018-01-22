/*
 *  TestTypes.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/19/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cassert>
#include <stdio.h>

#include "common.h"


struct StateVector
{
   Real a1r1, a2r2, ru, rv, rw, e, a2, dummy;

   void clear() {a1r1 = a2r2 = ru = rv = rw = e = a2 = dummy = 0.0;}
};

//union GP { StateVector s; StateVector dsdt; };
struct GP
{
	StateVector s; 
	StateVector dsdt;

	void clear()
	{
		s.clear();
		dsdt.clear();
	}
};

template < int _LX, int _LY, int _LZ> 
struct MatrixGP
{
	static const int LX = _LX;
	static const int LY = _LY;
	static const int LZ = _LZ;
	
	//char screwup_alignment;
	GP __attribute__((__aligned__(16))) data[_LZ][_LY][_LX];
	
	
	inline GP& operator()(const int ix, const int iy, const int iz)
	{
		assert(ix >= 0); assert(ix < _LX);
		assert(iy >= 0); assert(iy < _LY);
		assert(iy >= 0); assert(iy < _LY);
		
		return data[iz][iy][ix];
	}
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(MatrixGP)/1024.;
	}
	
	void compare(MatrixGP& block, double accuracy, string kernelname, bool compare_dsdt=true)
	{
		double maxe[7]={0,0,0,0,0,0,0};
		double sume[7]={0,0,0,0,0,0,0};
		
		for(int iz = 0; iz<_BLOCKSIZE_; iz++)
			for(int iy = 0; iy<_BLOCKSIZE_; iy++)
				for(int ix = 0; ix<_BLOCKSIZE_; ix++)
				{
					StateVector a = (*this)(ix, iy, iz).dsdt;
					StateVector b = block(ix, iy, iz).dsdt;
					
					//if we should not compare "dsdt" then we compare "s"
					if (!compare_dsdt)
					{
						a = (*this)(ix, iy, iz).s;
						b = block(ix, iy, iz).s;
					}
                    
                    const double s[7]  = {
						b.a1r1 ,
                                                b.a2r2,
						b.ru ,
						b.rv ,
						b.rw ,
						b.e ,
						b.a2
					};
					
					const double e[7]  = {
						b.a1r1 - a.a1r1,
						b.ru - a.ru,
						b.rv - a.rv,
						b.rw - a.rw,
						b.e - a.e,
						b.a2 - a.a2
					};
					
					for(int i=0; i<7; ++i)
						assert(!isnan(e[i]));
					
					for(int i=0; i<7; i++)
						if (fabs(e[i])/fabs(s[i])>accuracy && fabs(e[i])>accuracy) printf("significant error at %d %d %d %d -> e=%e (rel is %e, values are %e %e)\n", ix, iy, iz, i, e[i], e[i]/s[i],s[i],s[i]-e[i]);
					
					for(int i=0; i<7; i++)
						maxe[i] = max(fabs(e[i]), maxe[i]);
					
					for(int i=0; i<7; i++)
						sume[i] += fabs(e[i]);
                    
                    printf("\tLinf discrepancy:\t");
                    for(int i=0; i<7; i++)
                        printf("%.2e ", maxe[i]);
                    
                    printf("\n\tL1 (dh=1):       \t");
                    for(int i=0; i<7; i++)
                        printf("%.2e ", sume[i]);
                    printf("\n");
				}
	}
};

typedef MatrixGP<_BLOCKSIZE_, _BLOCKSIZE_, _BLOCKSIZE_> Block;

template<int WIDTH>
struct TestLab : MatrixGP<_BLOCKSIZE_+2*WIDTH, _BLOCKSIZE_+2*WIDTH, _BLOCKSIZE_+2*WIDTH>
{
        inline GP& operator()(const int ix, const int iy, const int iz)
        {
                assert(ix >= -WIDTH); assert(ix < _BLOCKSIZE_+WIDTH);
                assert(iy >= -WIDTH); assert(iy < _BLOCKSIZE_+WIDTH);
                assert(iz >= -WIDTH); assert(iz < _BLOCKSIZE_+WIDTH);
                
                return this->data[iz+WIDTH][iy+WIDTH][ix+WIDTH];
        }
        
        static float kB(int nobjects=1)
        {
                return nobjects*sizeof(TestLab)/1024.;
        }
};


/*
struct TestLab : MatrixGP<_BLOCKSIZE_+6, _BLOCKSIZE_+6, _BLOCKSIZE_+6>
{
	inline GP& operator()(const int ix, const int iy, const int iz)
	{
		assert(ix >= -3); assert(ix < _BLOCKSIZE_+3);
		assert(iy >= -3); assert(iy < _BLOCKSIZE_+3);
		assert(iz >= -3); assert(iz < _BLOCKSIZE_+3);
		
		return data[iz+3][iy+3][ix+3];
	}
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(TestLab)/1024.;
	}
};

struct TestLab_S2 : MatrixGP<_BLOCKSIZE_+2, _BLOCKSIZE_+2, _BLOCKSIZE_+2>
{
	inline GP& operator()(const int ix, const int iy, const int iz)
	{
		assert(ix >= -1); assert(ix < _BLOCKSIZE_+1);
		assert(iy >= -1); assert(iy < _BLOCKSIZE_+1);
		assert(iz >= -1); assert(iz < _BLOCKSIZE_+1);
		
		return data[iz+1][iy+1][ix+1];
	}
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(TestLab_S2)/1024.;
	}
};
*/

