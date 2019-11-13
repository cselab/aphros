/*
 *  GridMPI.h
 *  Cubism (allocation free)
 *
 *  Created by Fabian Wermelinger on 04/11/19.
 *  Copyright 2019 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#ifndef GRIDMPI_H_J17BPEDT
#define GRIDMPI_H_J17BPEDT

#include <map>
#include <mpi.h>
#include <string>
#include <vector>

#include "Cubism/StencilInfo.h"
#include "Cubism/afree/SynchronizerMPI.h"

template <typename TGrid>
class GridMPI : public TGrid
{
public:
    using Block = typename TGrid::Block;
    using BlockType = typename TGrid::BlockType;
    using Synch = SynchronizerMPI<Block, Block::bx, Block::by, Block::bz>;
    friend Synch;

private:
    size_t timestamp;

protected:
    int myrank, mypeindex[3], pesize[3];
    int periodic[3];
    int mybpd[3], myblockstotalsize, blocksize[3];

    std::vector<BlockInfo> cached_blockinfo;

    // acts like cache
    std::map<std::string, Synch *> SynchronizerMPIs;

    MPI_Comm worldcomm;
    MPI_Comm cartcomm;

public:
    GridMPI(const int npeX,
            const int npeY,
            const int npeZ,
            const int nX,
            const int nY = 1,
            const int nZ = 1,
            const double maxextent = 1,
            const MPI_Comm comm = MPI_COMM_WORLD)
        : TGrid(nX, nY, nZ, maxextent), timestamp(0), worldcomm(comm)
    {
        blocksize[0] = Block::sizeX;
        blocksize[1] = Block::sizeY;
        blocksize[2] = Block::sizeZ;

        mybpd[0] = nX;
        mybpd[1] = nY;
        mybpd[2] = nZ;
        myblockstotalsize = nX * nY * nZ;

        periodic[0] = true;
        periodic[1] = true;
        periodic[2] = true;

        pesize[0] = npeX;
        pesize[1] = npeY;
        pesize[2] = npeZ;

        int world_size;
        MPI_Comm_size(worldcomm, &world_size);
        assert(npeX * npeY * npeZ == world_size);

        MPI_Cart_create(worldcomm, 3, pesize, periodic, true, &cartcomm);
        MPI_Comm_rank(cartcomm, &myrank);
        MPI_Cart_coords(cartcomm, myrank, 3, mypeindex);

        const std::vector<BlockInfo> vInfo = TGrid::getBlocksInfo();

        // cells total
        const int nn[3] = {mybpd[0] * blocksize[0] * pesize[0],
                           mybpd[1] * blocksize[1] * pesize[1],
                           mybpd[2] * blocksize[2] * pesize[2]};
        // cell size (h_gridpoint from BlockInfo)
        const double hc = (maxextent / std::max(nn[0], std::max(nn[1], nn[2])));
        // block extent
        const double h =
            std::max(blocksize[0], std::max(blocksize[1], blocksize[2])) * hc;

        for (size_t i = 0; i < vInfo.size(); ++i) {
            BlockInfo info = vInfo[i];
            info.h = h;
            info.h_gridpoint = hc;
            for (int j = 0; j < 3; ++j) {
                info.index[j] += mypeindex[j] * mybpd[j];
                info.origin[j] = info.index[j] * h;
            }
            cached_blockinfo.push_back(info);
        }
    }

    ~GridMPI()
    {
        for (typename std::map<std::string, Synch *>::const_iterator it =
                 SynchronizerMPIs.begin();
             it != SynchronizerMPIs.end();
             ++it)
            delete it->second;

        SynchronizerMPIs.clear();
        MPI_Comm_free(&cartcomm);
    }

    std::vector<BlockInfo> &getBlocksInfo() override
    {
        return cached_blockinfo;
    }

    const std::vector<BlockInfo> &getBlocksInfo() const override
    {
        return cached_blockinfo;
    }

    std::vector<BlockInfo> &getResidentBlocksInfo()
    {
        return TGrid::getBlocksInfo();
    }

    const std::vector<BlockInfo> &getResidentBlocksInfo() const
    {
        return TGrid::getBlocksInfo();
    }

    // for a given kernel p (e.g. Processing=Diffusion)
    // based on StencilInfo (bounding box and selected components)
    // returns a SynchronizerMPI (new or existing)
    // and performs communication
    // (e.g. two kernels with identical stencils and
    // same elected components would share a common SynchronizerMPI)
    template <typename Processing, typename TView>
    Synch &sync(const std::string &name,
                Processing &p,
                std::vector<std::vector<TView>> &fields)
    {
        const StencilInfo stencil = p.stencil;
        assert(stencil.isvalid());

        Synch *queryresult = NULL;

        typename std::map<std::string, Synch *>::iterator itSynchronizerMPI =
            SynchronizerMPIs.find(name);

        if (itSynchronizerMPI == SynchronizerMPIs.end()) {
            queryresult = new Synch(fields,
                                    getBlocksInfo(),
                                    stencil,
                                    cartcomm,
                                    mybpd,
                                    blocksize);

            SynchronizerMPIs[name] = queryresult;
        } else
            queryresult = itSynchronizerMPI->second;

        // perform communication
        queryresult->sync(sizeof(Real) > 4 ? MPI_DOUBLE : MPI_FLOAT, timestamp);

        timestamp = (timestamp + 1) % 32768;

        return *queryresult;
    }

    inline bool isRegistered(const std::string &name) const
    {
        return SynchronizerMPIs.find(name) != SynchronizerMPIs.end();
    }

    Synch &getSynchronizerMPI(const std::string &name)
    {
        assert((SynchronizerMPIs.find(name) != SynchronizerMPIs.end()));

        return *SynchronizerMPIs.find(name)->second;
    }

    size_t getResidentBlocksPerDimension(int idim) const
    {
        assert(idim >= 0 && idim < 3);
        return (size_t)mybpd[idim];
    }

    size_t getBlocksPerDimension(int idim) const override
    {
        assert(idim >= 0 && idim < 3);
        return (size_t)mybpd[idim] * pesize[idim];
    }

    void peindex(int mypeindex[3]) const
    {
        for (int i = 0; i < 3; ++i)
            mypeindex[i] = this->mypeindex[i];
    }

    size_t getTimeStamp() const { return timestamp; }

    MPI_Comm getCartComm() const { return cartcomm; }

    MPI_Comm getWorldComm() const { return worldcomm; }
};

#endif /* GRIDMPI_H_J17BPEDT */
