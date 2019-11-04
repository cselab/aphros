/*
 *  Grid_new.h
 *  Cubism (does not own block data)
 *
 *  Created by Fabian Wermelinger on 04/11/19.
 *  Copyright 2019 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cassert>
#include <vector>

#include "BlockInfo.h"

template <typename B_>
class Grid
{
public:
    using Block = B_;
    using BlockType = Block;

private:
    // FIXME: [fabianw@mavt.ethz.ch; 2019-11-04] BlockInfo type could be
    // replaced with templated type.
    std::vector<BlockInfo> m_vInfo;

protected:
    const size_t NX, NY, NZ, N;
    const double maxextent;

    size_t _encode(const size_t ix, const size_t iy, const size_t iz) const
    {
        assert(ix >= 0 && ix < NX);
        assert(iy >= 0 && iy < NY);
        assert(iz >= 0 && iz < NZ);

        return ix + NX * (iy + NY * iz);
    }

public:
    Grid(const size_t NX_,
         const size_t NY_ = 1,
         const size_t NZ_ = 1,
         const double maxextent_ = 1)
        : NX(NX_), NY(NY_), NZ(NZ_), N(NX_ * NY_ * NZ_), maxextent(maxextent_)
    {
        // FIXME: [fabianw@mavt.ethz.ch; 2019-11-04] Get a const reference to
        // external vector of block infos.  Requires that _encode() id is
        // correct in vector
        m_vInfo.reserve(N);
        const double o[3] = {0, 0, 0};
        for (size_t iz = 0; iz < NZ; iz++)
            for (size_t iy = 0; iy < NY; iy++)
                for (size_t ix = 0; ix < NX; ix++) {
                    const long long blockID = _encode(ix, iy, iz);
                    const int idx[3] = {(int)ix, (int)iy, (int)iz};
                    m_vInfo.push_back(BlockInfo(blockID, idx, o, 0, 0));
                }
    }

    virtual ~Grid() {}

    virtual size_t getBlocksPerDimension(int idim) const
    {
        assert(idim >= 0 && idim < 3);

        switch (idim) {
        case 0:
            return NX;
        case 1:
            return NY;
        case 2:
            return NZ;
        default:
            abort();
            return 0;
        }
    }

    virtual bool avail(int, int, int) const { return true; }

    virtual size_t operator()(size_t ix, size_t iy = 0, size_t iz = 0) const
    {
        return _encode((ix + NX) % NX, (iy + NY) % NY, (iz + NZ) % NZ);
    }

    virtual std::vector<BlockInfo> &getBlocksInfo() { return m_vInfo; }

    virtual const std::vector<BlockInfo> &getBlocksInfo() const
    {
        return m_vInfo;
    }
};
