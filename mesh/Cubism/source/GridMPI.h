/*
 *  GridMPI.h
 *  FacesMPI
 *
 *  Created by Diego Rossinelli on 10/21/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <vector>
#include <map>
#include <mpi.h>

#include "BlockInfo.h"
#include "StencilInfo.h"
#include "SynchronizerMPI.h"

template < typename TGrid >
class GridMPI : public TGrid
{
	size_t timestamp;

protected:

	friend class SynchronizerMPI;

	int myrank, mypeindex[3], pesize[3];
	int periodic[3];
	int mybpd[3], myblockstotalsize, blocksize[3];

	std::vector<BlockInfo> cached_blockinfo;

    std::map<StencilInfo, SynchronizerMPI *> SynchronizerMPIs;

    MPI_Comm worldcomm;
	MPI_Comm cartcomm;

public:

    typedef typename TGrid::BlockType Block;

	GridMPI(const int npeX, const int npeY, const int npeZ,
			const int nX, const int nY=1, const int nZ=1,
			const double maxextent = 1, const MPI_Comm comm = MPI_COMM_WORLD): TGrid(nX, nY, nZ, maxextent), timestamp(0), worldcomm(comm)
	{
		blocksize[0] = Block::sizeX;
		blocksize[1] = Block::sizeY;
		blocksize[2] = Block::sizeZ;

		mybpd[0] = nX;
		mybpd[1] = nY;
		mybpd[2] = nZ;
		myblockstotalsize = nX*nY*nZ;

		periodic[0] = true;
		periodic[1] = true;
		periodic[2] = true;

		pesize[0] = npeX;
		pesize[1] = npeY;
		pesize[2] = npeZ;

        int world_size;
        MPI_Comm_size(worldcomm, &world_size);
		assert(npeX*npeY*npeZ == world_size);

        MPI_Cart_create(worldcomm, 3, pesize, periodic, true, &cartcomm);
        MPI_Comm_rank(cartcomm, &myrank);
        MPI_Cart_coords(cartcomm, myrank, 3, mypeindex);

        for (int i = 0; i < 3; ++i)
            delete this->m_mesh_maps[i];
        std::vector<MeshMap<Block>*> tmp;
        this->m_mesh_maps.swap(tmp);

        const double blockwidth = (maxextent / std::max(nX*npeX, std::max(nY*npeY, nZ*npeZ)));
        const double extents[3] = {blockwidth*nX*npeX, blockwidth*nY*npeY, blockwidth*nZ*npeZ};
        const unsigned int nBlocks[3] = {nX*npeX, nY*npeY, nZ*npeZ};
        for (int i = 0; i < 3; ++i)
        {
            MeshMap<Block>* m = new MeshMap<Block>(0.0, extents[i], nBlocks[i]);
            UniformDensity uniform;
            m->init(&uniform); // uniform only for this constructor
            this->m_mesh_maps.push_back(m);
        }

		const std::vector<BlockInfo> vInfo = TGrid::getBlocksInfo();

        const double h_gridpoint = maxextent / (double)max(getBlocksPerDimension(0)*blocksize[0], max(getBlocksPerDimension(1)*blocksize[1], getBlocksPerDimension(2)*blocksize[2]));

        for(int i=0; i<vInfo.size(); ++i)
        {
            BlockInfo info = vInfo[i];

			info.h_gridpoint = h_gridpoint;

            info.h = info.h_gridpoint * blocksize[0];// only for blocksize[0]=blocksize[1]=blocksize[2]

            for(int j=0; j<3; ++j)
            {
                info.index[j] += mypeindex[j]*mybpd[j];

                info.origin[j] = this->m_mesh_maps[j]->block_origin(info.index[j]);

                info.uniform_grid_spacing[j] = h_gridpoint;

                info.block_extent[j] = this->m_mesh_maps[j]->block_width(info.index[j]);

                info.ptr_grid_spacing[j] = this->m_mesh_maps[j]->get_grid_spacing(info.index[j]);
            }

            cached_blockinfo.push_back(info);
        }
	}


    GridMPI(MeshMap<Block>* const mapX, MeshMap<Block>* const mapY, MeshMap<Block>* const mapZ,
            const int npeX, const int npeY, const int npeZ,
            const int nX=0, const int nY=0, const int nZ=0,
            const MPI_Comm comm = MPI_COMM_WORLD): TGrid(mapX,mapY,mapZ,nX,nY,nZ), timestamp(0), worldcomm(comm)
    {
        blocksize[0] = Block::sizeX;
        blocksize[1] = Block::sizeY;
        blocksize[2] = Block::sizeZ;

        assert(mapX->nblocks() == npeX*nX);
        assert(mapY->nblocks() == npeY*nY);
        assert(mapZ->nblocks() == npeZ*nZ);
        mybpd[0] = nX;
        mybpd[1] = nY;
        mybpd[2] = nZ;
        myblockstotalsize = nX*nY*nZ;

        periodic[0] = true;
        periodic[1] = true;
        periodic[2] = true;

        pesize[0] = npeX;
        pesize[1] = npeY;
        pesize[2] = npeZ;

        int world_size;
        MPI_Comm_size(worldcomm, &world_size);
        assert(npeX*npeY*npeZ == world_size);

        MPI_Cart_create(worldcomm, 3, pesize, periodic, true, &cartcomm);
        MPI_Comm_rank(cartcomm, &myrank);
        MPI_Cart_coords(cartcomm, myrank, 3, mypeindex);

        const std::vector<BlockInfo> vInfo = TGrid::getBlocksInfo();

        MeshMap<Block>* const ptr_map[3] = {mapX, mapY, mapZ};
        for(int i=0; i<vInfo.size(); ++i)
        {
            BlockInfo info = vInfo[i];

            for(int j=0; j<3; ++j)
            {
                info.index[j] += mypeindex[j]*mybpd[j];

                info.origin[j] = ptr_map[j]->block_origin(info.index[j]);

                info.block_extent[j] = ptr_map[j]->block_width(info.index[j]);

                info.ptr_grid_spacing[j] = ptr_map[j]->get_grid_spacing(info.index[j]);
            }

            cached_blockinfo.push_back(info);
        }
    }


	~GridMPI()
	{
		for( std::map<StencilInfo, SynchronizerMPI*>::const_iterator it = SynchronizerMPIs.begin(); it != SynchronizerMPIs.end(); ++it)
			delete it->second;

		SynchronizerMPIs.clear();
        MPI_Comm_free(&cartcomm);
	}

	std::vector<BlockInfo>& getBlocksInfo()
	{
		return cached_blockinfo;
	}

	const std::vector<BlockInfo>& getBlocksInfo() const
	{
		return cached_blockinfo;
	}

	std::vector<BlockInfo>& getResidentBlocksInfo()
	{
		return TGrid::getBlocksInfo();
	}

	const std::vector<BlockInfo>& getResidentBlocksInfo() const
	{
		return TGrid::getBlocksInfo();
	}

	virtual bool avail(int ix, int iy=0, int iz=0) const
	{
		//return true;
		const int originX = mypeindex[0]*mybpd[0];
		const int originY = mypeindex[1]*mybpd[1];
		const int originZ = mypeindex[2]*mybpd[2];

		const int nX = pesize[0]*mybpd[0];
		const int nY = pesize[1]*mybpd[1];
		const int nZ = pesize[2]*mybpd[2];

		ix = (ix + nX) % nX;
		iy = (iy + nY) % nY;
		iz = (iz + nZ) % nZ;

		const bool xinside = (ix>= originX && ix<nX);
		const bool yinside = (iy>= originY && iy<nY);
		const bool zinside = (iz>= originZ && iz<nZ);

		assert(TGrid::avail(ix-originX, iy-originY, iz-originZ));
		return xinside && yinside && zinside;
	}

	inline Block& operator()(int ix, int iy=0, int iz=0) const
	{
		//assuming ix,iy,iz to be global
		const int originX = mypeindex[0]*mybpd[0];
		const int originY = mypeindex[1]*mybpd[1];
		const int originZ = mypeindex[2]*mybpd[2];

		const int nX = pesize[0]*mybpd[0];
		const int nY = pesize[1]*mybpd[1];
		const int nZ = pesize[2]*mybpd[2];

		ix = (ix + nX) % nX;
		iy = (iy + nY) % nY;
		iz = (iz + nZ) % nZ;

		assert(ix>= originX && ix<nX);
		assert(iy>= originY && iy<nY);
		assert(iz>= originZ && iz<nZ);

		return TGrid::operator()(ix-originX, iy-originY, iz-originZ);
	}

	template<typename Processing>
	SynchronizerMPI& sync(Processing& p)
	{
		const StencilInfo stencil = p.stencil;
		assert(stencil.isvalid());

		SynchronizerMPI * queryresult = NULL;

		typename std::map<StencilInfo, SynchronizerMPI*>::iterator itSynchronizerMPI = SynchronizerMPIs.find(stencil);

		if (itSynchronizerMPI == SynchronizerMPIs.end())
		{
			queryresult = new SynchronizerMPI(SynchronizerMPIs.size(), stencil, getBlocksInfo(), cartcomm, mybpd, blocksize);

			SynchronizerMPIs[stencil] = queryresult;
		}
		else  queryresult = itSynchronizerMPI->second;

		queryresult->sync(sizeof(typename Block::element_type)/sizeof(Real), sizeof(Real)>4 ? MPI_DOUBLE : MPI_FLOAT, timestamp);

		timestamp = (timestamp + 1) % 32768;

		return *queryresult;
	}

	template<typename Processing>
	const SynchronizerMPI& get_SynchronizerMPI(Processing& p) const
	{
		assert((SynchronizerMPIs.find(p.stencil) != SynchronizerMPIs.end()));

		return *SynchronizerMPIs.find(p.stencil)->second;
	}

	int getResidentBlocksPerDimension(int idim) const
	{
		assert(idim>=0 && idim<3);
		return mybpd[idim];
	}

	int getBlocksPerDimension(int idim) const
	{
		assert(idim>=0 && idim<3);
		return mybpd[idim]*pesize[idim];
	}

	void peindex(int mypeindex[3]) const
	{
		for(int i=0; i<3; ++i)
			mypeindex[i] = this->mypeindex[i];
	}

    size_t getTimeStamp() const
    {
        return timestamp;
    }

    MPI_Comm getCartComm() const
	{
		return cartcomm;
	}

    MPI_Comm getWorldComm() const
	{
		return worldcomm;
	}
};
