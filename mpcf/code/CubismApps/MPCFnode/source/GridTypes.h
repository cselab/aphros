/*
 *  GridTypes.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger on 05/10/17
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef GRIDTYPES_H_LXM9TWTK
#define GRIDTYPES_H_LXM9TWTK

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif


struct FluidElement
{
  Real alpha2;
  Real a[7];

    void clear() { alpha2 = 0.0; }

    void init(const Real val) { alpha2 = val; }

    FluidElement& operator = (const FluidElement & gp)
    {
        this->alpha2 = gp.alpha2;
        return *this;
    }
};

FluidElement operator * (const Real a, FluidElement gp);
FluidElement operator + (FluidElement gpa, FluidElement gpb);
FluidElement operator - (FluidElement gpa, FluidElement gpb);


struct FluidBlock
{
    static const int sizeX;
    static const int sizeY;
    static const int sizeZ;

    static const int gptfloats = sizeof(FluidElement)/sizeof(Real);

    typedef FluidElement ElementType;
    typedef FluidElement element_type;

    FluidElement __attribute__((__aligned__(_ALIGNBYTES_))) data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];

    Real __attribute__((__aligned__(_ALIGNBYTES_))) tmp[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][gptfloats];

    void clear_data()
    {
        const int N = sizeX*sizeY*sizeZ;
        FluidElement * const e = &data[0][0][0];
        for(int i=0; i<N; ++i) e[i].clear();
    }

    void clear_tmp()
    {
        const int N = sizeX * sizeY * sizeZ * gptfloats;

        Real * const e = &tmp[0][0][0][0];
        for(int i=0; i<N; ++i) e[i] = 0;
    }

    void clear()
    {
        clear_data();
        clear_tmp();
    }

    inline FluidElement& operator()(int ix, int iy=0, int iz=0)
    {
        assert(ix>=0 && ix<sizeX);
        assert(iy>=0 && iy<sizeY);
        assert(iz>=0 && iz<sizeZ);

        return data[iz][iy][ix];
    }

    template <typename Streamer>
    inline void Write(ofstream& output, Streamer streamer) const
    {
        for(int iz=0; iz<sizeZ; iz++)
            for(int iy=0; iy<sizeY; iy++)
                for(int ix=0; ix<sizeX; ix++)
                    streamer.operate(data[iz][iy][ix], output);
    }

    template <typename Streamer>
    inline void Read(ifstream& input, Streamer streamer)
    {
        for(int iz=0; iz<sizeZ; iz++)
            for(int iy=0; iy<sizeY; iy++)
                for(int ix=0; ix<sizeX; ix++)
                    streamer.operate(input, data[iz][iy][ix]);
    }
};


typedef FluidBlock Block_t;  
typedef Grid<Block_t, std::allocator> GridBase;
typedef GridBase Grid_t;

#endif /* GRIDTYPES_H_LXM9TWTK */
