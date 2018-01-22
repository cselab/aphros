

/*
 *  *  *  SOA2D.h
 *   *   *
 *    *    *
 *     *     *  Created by Diego Rossinelli on 5/15/12.
 *      *      *  Copyright 2012 ETH Zurich. All rights reserved.
 *       *       *
 *        *        */

#pragma once

#include <cassert>

#include "common.h"

template < int _SX, int _EX, int _SY, int _EY, typename TReal=Real >
#ifdef __xlC__
__align(_ALIGNBYTES_)
#endif
struct SOA2D
{
	static const int _CPERALIGNBYTES = _ALIGNBYTES_/sizeof(TReal);

	static const int SX = _CPERALIGNBYTES*((_SX - (_CPERALIGNBYTES-1))/_CPERALIGNBYTES);
	static const int EX = _CPERALIGNBYTES*((_EX + (_CPERALIGNBYTES-1))/_CPERALIGNBYTES);
	static const int NX = _EX - _SX;
	static const int NY = _EY - _SY;

	static const int PITCH = EX - SX;

	__attribute__((aligned(_ALIGNBYTES_))) TReal data[NY][PITCH];


	SOA2D()
	{
		assert(((size_t)(&data[0][0]) % _ALIGNBYTES_) == 0);
	}

	inline TReal operator()(const int ix, const int iy) const
	{
		assert(ix >= SX); assert(ix < EX);
		assert(iy >= _SY); assert(iy < _EY);

		return data[iy-_SY][ix-SX];
	}

	inline const TReal * ptr(const int ix, const int iy) const
	{
		assert(ix >= SX); assert(ix < EX);
		assert(iy >= _SY); assert(iy < _EY);

		return &data[iy-_SY][ix-SX];
	}

	inline TReal& ref(const int ix, const int iy)
	{
		assert(ix >= SX); assert(ix < EX);
		assert(iy >= _SY); assert(iy < _EY);

		return data[iy-_SY][ix-SX];
	}

	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(SOA2D)/1024.;
	}

	inline SOA2D& operator= (const SOA2D& c)
	{
		for(int y=0; y<NY; ++y)
			for(int x=0; x<PITCH; ++x)
				data[y][x] = c.data[y][x];

		return *this;
	}

        inline SOA2D& operator-= (const SOA2D& c)
        {
                for(int y=0; y<NY; ++y)
                        for(int x=0; x<PITCH; ++x)
                                data[y][x] -= c.data[y][x];

                return *this;
        }

        inline SOA2D& operator/= (const SOA2D& c)
        {
                for(int y=0; y<NY; ++y)
                        for(int x=0; x<PITCH; ++x)
                                data[y][x] /= c.data[y][x];

                return *this;
        }

        inline SOA2D& operator+= (const SOA2D& c)
        {
                for(int y=0; y<NY; ++y)
                        for(int x=0; x<PITCH; ++x)
                                data[y][x] += c.data[y][x];

                return *this;
        }

        inline SOA2D& operator*= (const SOA2D& c)
        {
                for(int y=0; y<NY; ++y)
                        for(int x=0; x<PITCH; ++x)
                                data[y][x] *= c.data[y][x];

                return *this;
        }

        inline SOA2D& operator+= (const Real c)
        {
                for(int y=0; y<NY; ++y)
                        for(int x=0; x<PITCH; ++x)
                                data[y][x] += c;

                return *this;
        }

        inline SOA2D& operator-= (const Real c)
        {
                for(int y=0; y<NY; ++y)
                        for(int x=0; x<PITCH; ++x)
                                data[y][x] -= c;

                return *this;
        }

        inline SOA2D& operator*= (const Real c)
        {
                for(int y=0; y<NY; ++y)
                        for(int x=0; x<PITCH; ++x)
                                data[y][x] *= c;

                return *this;
        }

};

template< int _SX, int _EX, int _SY, int _EY, int _NSLICES, typename TReal=Real>
struct RingSOA2D
{
	int currslice;

	SOA2D<_SX, _EX, _SY, _EY, TReal> slices[_NSLICES];

	RingSOA2D(): currslice(0){ }

	inline const SOA2D<_SX, _EX, _SY, _EY, TReal>& operator()(const int relativeid=0) const
	{
		return slices[(relativeid + currslice + _NSLICES) % _NSLICES];
	}

	inline SOA2D<_SX, _EX, _SY, _EY, TReal>& ref(const int relativeid=0)
	{
		return slices[(relativeid + currslice + _NSLICES) % _NSLICES];
	}

	void next(){ currslice = (currslice + 1) % _NSLICES; }

	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(RingSOA2D)/1024.;
	}
};

typedef SOA2D<0, _BLOCKSIZE_, 0, _BLOCKSIZE_> OutputSOA;

/*working dataset types */
typedef SOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3> InputSOA;
typedef RingSOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3, 6> RingInputSOA;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_> TempSOA;
typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_, 2> RingTempSOA;
typedef RingSOA2D<0,_BLOCKSIZE_+1, 0, _BLOCKSIZE_, 3> RingTempSOA3;
