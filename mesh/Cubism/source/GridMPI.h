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

using namespace std;

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
	bool periodic[3];
	int mybpd[3], myblockstotalsize, blocksize[3];

	vector<BlockInfo> cached_blockinfo;

	map<StencilInfo, SynchronizerMPI *> SynchronizerMPIs;

	MPI::Cartcomm cartcomm;

public:

    typedef typename TGrid::BlockType Block;

	GridMPI(const int npeX, const int npeY, const int npeZ,
			const int nX, const int nY=1, const int nZ=1,
			const double maxextent = 1): TGrid(nX, nY, nZ, maxextent), timestamp(0)
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

		assert(npeX*npeY*npeZ == MPI::COMM_WORLD.Get_size());

		cartcomm = MPI::COMM_WORLD.Create_cart(3, pesize, periodic, true);
		myrank = cartcomm.Get_rank();

		cartcomm.Get_coords(myrank, 3, mypeindex);

		const vector<BlockInfo> vInfo = TGrid::getBlocksInfo();

		for(int i=0; i<vInfo.size(); ++i)
		{
			BlockInfo info = vInfo[i];

			info.h_gridpoint = maxextent / (double)max (getBlocksPerDimension(0)*blocksize[0],
														max(getBlocksPerDimension(1)*blocksize[1],
															getBlocksPerDimension(2)*blocksize[2]));

            info.h = info.h_gridpoint * blocksize[0];// only for blocksize[0]=blocksize[1]=blocksize[2]

            for(int j=0; j<3; ++j)
			{
				info.index[j] += mypeindex[j]*mybpd[j];
				info.origin[j] = info.index[j]*info.h;
			}

			cached_blockinfo.push_back(info);
		}
	}

	~GridMPI()
	{
		for( map<StencilInfo, SynchronizerMPI*>::const_iterator it = SynchronizerMPIs.begin(); it != SynchronizerMPIs.end(); ++it)
			delete it->second;

		SynchronizerMPIs.clear();
	}

	vector<BlockInfo>& getBlocksInfo()
	{
		return cached_blockinfo;
	}

	const vector<BlockInfo>& getBlocksInfo() const
	{
		return cached_blockinfo;
	}

	vector<BlockInfo>& getResidentBlocksInfo()
	{
		return TGrid::getBlocksInfo();
	}

	const vector<BlockInfo>& getResidentBlocksInfo() const
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

		typename map<StencilInfo, SynchronizerMPI*>::iterator itSynchronizerMPI = SynchronizerMPIs.find(stencil);

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

    MPI::Cartcomm getCartComm() const
	{
		return cartcomm;
	}

    double getH() const
    {
        vector<BlockInfo> vInfo = this->getBlocksInfo();
        BlockInfo info = vInfo[0];
        return info.h_gridpoint;
    }
};
