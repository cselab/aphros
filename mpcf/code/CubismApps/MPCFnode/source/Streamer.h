/*
 *  Streamer.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 09/13/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef STREAMER_H_0HQKKPXE
#define STREAMER_H_0HQKKPXE

#include <cassert>
#include <cmath>
#include <string>
#include <algorithm>
#include "Types.h"

struct StreamerAlpha2
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    const Block_t& ref;

    StreamerAlpha2(const Block_t& b): ref(b) {}

    inline void operate(const int ix, const int iy, const int iz, Real output[NCHANNELS]) const
    {
        const FluidElement& el = ref.data[iz][iy][ix];
        output[0] = el.alpha2;
    }

    inline Real operate(const int ix, const int iy, const int iz) const
    {
        const FluidElement& el = ref.data[iz][iy][ix];
        return el.alpha2;
    }

    static const char * getAttributeName() { return "Scalar"; }
};


#endif /* STREAMER_H_0HQKKPXE */
