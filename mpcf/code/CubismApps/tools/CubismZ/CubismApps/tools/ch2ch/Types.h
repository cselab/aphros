/*
 *  Types.h
 *  MPCFnode
 *
 *  Created by Jonas Sukys on May 10, 2015.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include <fstream>
#include "math.h"

using namespace std;

#include <Grid.h>

class Simulation
{
public:
    
	virtual void setup() { }
	virtual void dispose() { }
	virtual ~Simulation() { }
};

struct FluidElement
{
  Real r, u, v, w, p, G, P, dummy;
  
  void clear() { r = u = v = w = p = G = P = dummy = 0; }
  
  FluidElement& operator = (const FluidElement & gp)
  {
    this->r = gp.r;
    this->u = gp.u;
    this->v = gp.v;
    this->w = gp.w;
    this->p = gp.p;
    this->G = gp.G;
    this->P = gp.P;
    
    return *this;
  }
};

struct FluidBlock
{
	static const int sizeX = _BLOCKSIZE_;
	static const int sizeY = _BLOCKSIZE_;
	static const int sizeZ = _BLOCKSIZE_;

	static const int gptfloats = sizeof(FluidElement)/sizeof(Real);
	
	typedef FluidElement ElementType;
	typedef FluidElement element_type;
	
	FluidElement __attribute__((__aligned__(_ALIGNBYTES_))) data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
    
	void clear_data()
	{
		const int N = sizeX*sizeY*sizeZ;
		FluidElement * const e = &data[0][0][0];
		for(int i=0; i<N; ++i) e[i].clear();
	}
    
	void clear()
	{
		clear_data();
	}
    
	inline FluidElement& operator()(int ix, int iy=0, int iz=0)
	{
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
		assert(iz>=0 && iz<sizeZ);
		
		return data[iz][iy][ix];
	}
		
	template <typename Streamer>
	inline void minmax(Real minval[Streamer::channels], Real maxval[Streamer::channels], Streamer streamer = Streamer())
	{
		enum { NCHANNELS = Streamer::channels };
				
		streamer.operate(data[0][0][0], minval);
		streamer.operate(data[0][0][0], maxval);
		
		for(int iz=0; iz<sizeZ; iz++)
			for(int iy=0; iy<sizeY; iy++)
				for(int ix=0; ix<sizeX; ix++)
				{
					Real tmp[NCHANNELS];
					
					streamer.operate(data[iz][iy][ix], tmp);
					
					for(int ic = 0; ic < NCHANNELS; ++ic)
						minval[ic] = std::min(minval[ic], tmp[ic]);
					
					for(int ic = 0; ic < NCHANNELS; ++ic)
						maxval[ic] = std::max(maxval[ic], tmp[ic]);
				}
	}
};

struct StreamerGridPointIterative
{
	static const int channels = 7;
  
	FluidBlock * ref;
	StreamerGridPointIterative(FluidBlock& b): ref(&b) {}
	StreamerGridPointIterative(): ref(NULL) {}

	template<int channel>
	static inline Real operate(const FluidElement& input) { abort(); return 0; }

	inline Real operate(const int ix, const int iy, const int iz) const
	{
		cout << "You must not call this operate method of StreamerGridPointIterative" << endl;
		abort();
		return 0;
	}

	const char * name() { return "StreamerGridPointIterative" ; }
};

template<> inline Real StreamerGridPointIterative::operate<0>(const FluidElement& e) { return e.r; }
template<> inline Real StreamerGridPointIterative::operate<1>(const FluidElement& e) { return e.u; }
template<> inline Real StreamerGridPointIterative::operate<2>(const FluidElement& e) { return e.v; }
template<> inline Real StreamerGridPointIterative::operate<3>(const FluidElement& e) { return e.w; }
template<> inline Real StreamerGridPointIterative::operate<4>(const FluidElement& e) { return e.p; }
template<> inline Real StreamerGridPointIterative::operate<5>(const FluidElement& e) { return e.G; }
template<> inline Real StreamerGridPointIterative::operate<6>(const FluidElement& e) { return e.P; }
template<> inline Real StreamerGridPointIterative::operate<7>(const FluidElement& e) { return e.dummy; }

struct StreamerVelocityMagnitude
{
  static const int channels = 1;
  
  FluidBlock * ref;
  StreamerVelocityMagnitude (FluidBlock& b) : ref (&b) {}
  StreamerVelocityMagnitude () : ref (NULL) {}
  
  template<int channel>
  static inline Real operate(const FluidElement& input) { abort(); return 0; }
  
  inline Real operate(const int ix, const int iy, const int iz) const
  {
    const Real& u = ref->data[iz][iy][ix].u;
    const Real& v = ref->data[iz][iy][ix].v;
    const Real& w = ref->data[iz][iy][ix].w;
    
    Real m = sqrt ( u*u + v*v + w*w );
    
    return m;
  }
  
  const char * name() { return "StreamerVelocityMagnitude" ; }
};

struct StreamerSpeedOfSound
{
  static const int channels = 1;
  
  FluidBlock * ref;
  StreamerSpeedOfSound (FluidBlock& b) : ref (&b) {}
  StreamerSpeedOfSound () : ref (NULL) {}
  
  template<int channel>
  static inline Real operate(const FluidElement& input) { abort(); return 0; }
  
  inline Real operate(const int ix, const int iy, const int iz) const
  {
    const Real& r = ref->data[iz][iy][ix].r;
    const Real& p = ref->data[iz][iy][ix].p;
    const Real& G = ref->data[iz][iy][ix].G;
    const Real& P = ref->data[iz][iy][ix].P;
    
    Real c = sqrt ( (1/G + 1) * (p + P/G / (1/G + 1)) / r);
    
    return c;
  }
  
  const char * name() { return "StreamerSpeedOfSound" ; }
};

struct StreamerMachNumber
{
  static const int channels = 1;
  
  FluidBlock * ref;
  StreamerMachNumber (FluidBlock& b) : ref (&b) {}
  StreamerMachNumber () : ref (NULL) {}
  
  template<int channel>
  static inline Real operate(const FluidElement& input) { abort(); return 0; }
  
  inline Real operate(const int ix, const int iy, const int iz) const
  {
    const Real& r = ref->data[iz][iy][ix].r;
    const Real& u = ref->data[iz][iy][ix].u;
    const Real& v = ref->data[iz][iy][ix].v;
    const Real& w = ref->data[iz][iy][ix].w;
    const Real& p = ref->data[iz][iy][ix].p;
    const Real& G = ref->data[iz][iy][ix].G;
    const Real& P = ref->data[iz][iy][ix].P;
    
    Real m = sqrt ( u*u + v*v + w*w );
    Real c = sqrt ( (1/G + 1) * (p + P/G / (1/G + 1)) / r);
    Real M = m / c;
    
    return M;
  }
  
  const char * name() { return "StreamerMachNumber" ; }
};


#if 1
struct StreamerDummy_HDF5raw
{
	static const int NCHANNELS = 7;

	FluidBlock& ref;

	StreamerDummy_HDF5raw(FluidBlock& b): ref(b){}

	void operate(const int ix, const int iy, const int iz, Real output[7]) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];

		output[0] = input.r;
		output[1] = input.u;
		output[2] = input.v;
		output[3] = input.w;
		output[4] = input.p;
		output[5] = input.G;
		output[6] = input.P;
	}

	void operate(const Real output[7], const int ix, const int iy, const int iz) const
	{
		FluidElement& input = ref.data[iz][iy][ix];

		input.r = output[0];
		input.u = output[1];
		input.v = output[2];
		input.w = output[3];
		input.p = output[4];
		input.G = output[5];
		input.P = output[6];
	}

	void operate(const int ix, const int iy, const int iz, Real *ovalue, const int field) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];

		switch(field) {
		case 0: *ovalue = input.r; break;
		case 1: *ovalue = input.u; break;
		case 2: *ovalue = input.v; break;
		case 3: *ovalue = input.w; break;
		case 4: *ovalue = input.p; break;
		case 5: *ovalue = input.G; break;
		case 6: *ovalue = input.P; break;
		default: printf("unknown field\n"); abort(); break;
		}
	}

	void operate(const Real ivalue, const int ix, const int iy, const int iz, const int field) const
	{
		FluidElement& input = ref.data[iz][iy][ix];

		switch(field) {
		case 0:  input.r       = ivalue; break;
		case 1:  input.u       = ivalue; break;
		case 2:  input.v       = ivalue; break;
		case 3:  input.w       = ivalue; break;
		case 4:  input.p       = ivalue; break;
		case 5:  input.G       = ivalue; break;
		case 6:  input.P       = ivalue; break;
		default: printf("unknown field\n"); abort(); break;
		}
	}

	static const char * getAttributeName() { return "Tensor"; }
};
#endif

typedef Grid<FluidBlock, std::allocator> FluidGridBase;
typedef FluidGridBase FluidGrid;


