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
  Real alpha1rho1, alpha2rho2, ru, rv, rw, energy, alpha2, dummy;

    void clear() { alpha1rho1 = alpha2rho2 = ru = rv = rw = energy = alpha2 = dummy = 0.0; }

    void init(const Real val) { alpha1rho1 = alpha2rho2 = ru = rv = rw = energy = alpha2 = dummy = val; }

    FluidElement& operator = (const FluidElement & gp)
    {
        this->alpha1rho1 = gp.alpha1rho1;
        this->alpha2rho2 = gp.alpha2rho2;
        this->ru = gp.ru;
        this->rv = gp.rv;
        this->rw = gp.rw;
        this->energy=gp.energy;
        this->alpha2 = gp.alpha2;
        this->dummy = gp.dummy;

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

// TODO: VOLFRAC_EQSYS_URSULA: include again
/*
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
*/
};


struct FluidBlockNonUniform
{
    static const int sizeX;
    static const int sizeY;
    static const int sizeZ;

    static const int gptfloats = sizeof(FluidElement)/sizeof(Real);

    typedef FluidElement ElementType;
    typedef FluidElement element_type;

    FluidElement __attribute__((__aligned__(_ALIGNBYTES_))) data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];

    Real __attribute__((__aligned__(_ALIGNBYTES_))) tmp[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][gptfloats];

    Coefficients_t __attribute__((__aligned__(_ALIGNBYTES_))) coeffs_x[2]; // minus/plus
    Coefficients_t __attribute__((__aligned__(_ALIGNBYTES_))) coeffs_y[2]; // minus/plus
    Coefficients_t __attribute__((__aligned__(_ALIGNBYTES_))) coeffs_z[2]; // minus/plus
    Real __attribute__((__aligned__(_ALIGNBYTES_))) invh_x[_BLOCKSIZE_];
    Real __attribute__((__aligned__(_ALIGNBYTES_))) invh_y[_BLOCKSIZE_];
    Real __attribute__((__aligned__(_ALIGNBYTES_))) invh_z[_BLOCKSIZE_];

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
};


// TODO: [fabianw@mavt.ethz.ch; Wed May 10 2017 06:58:52 PM (-0700)] This is
// similar to the BlockLab issue where we need a typedef for the Lab type.
// Currently we must live with this.
#ifdef _NONUNIFORM_BLOCK_
typedef FluidBlockNonUniform Block_t; // nonuniform block
#else
typedef FluidBlock Block_t;  // uniform block
#endif /* _NONUNIFORM_BLOCK_ */

typedef Grid<Block_t, std::allocator> GridBase;
typedef GridBase Grid_t;

#endif /* GRIDTYPES_H_LXM9TWTK */
