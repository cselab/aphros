/*
 *  GridMorton.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/24/11.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

using namespace std;

#include "Grid.h"
#include "Indexers.h"

template <typename TGrid>
class GridMorton: public TGrid
{
protected:

	vector<int> f2m, m2f;

	vector<BlockInfo> cached_infos;

	void _generate_morton_mapping()
	{
		const int N = this->N;
		const int NX = this->NX;
		const int NY = this->NY;
		const int NZ = this->NZ;
		const int MAXND = max(NX, max(NY, NZ));

		IndexerMorton indexer(MAXND, MAXND, MAXND);

		vector< pair< int, int > > tobesorted(N);

		for(int iflat = 0; iflat < N; ++iflat)
		{
			const int ix = iflat % NX;
			const int iy = (iflat / NX) % NY;
			const int iz = (iflat / (NX * NY)) % NZ;

			tobesorted[iflat] = make_pair(indexer.encode(ix, iy, iz), iflat);
		}

		std::sort(tobesorted.begin(), tobesorted.end());

		m2f.resize(N);
		for(int imorton = 0; imorton < N; ++imorton)
			m2f[imorton] = tobesorted[imorton].second;

		f2m.resize(N);
		for(int imorton = 0; imorton < N; ++imorton)
			f2m[tobesorted[imorton].second] = imorton;
	}

	vector<BlockInfo> _getBlocksInfo() const
	{
		const int N = this->N;
		const unsigned int NX = this->NX;
		const unsigned int NY = this->NY;
		const unsigned int NZ = this->NZ;

		std::vector<BlockInfo> r(N);

		const double h = (this->maxextent / max(NX, max(NY, NZ)));

		for(int imorton = 0; imorton < N; ++imorton)
		{
			const int iflat = m2f[imorton];

			const int ix = iflat % NX;
			const int iy = (iflat / NX) % NY;
			const int iz = (iflat / (NX * NY)) % NZ;

			const int idx[3] = {ix, iy, iz};
			const double origin[3] = { ix * h, iy * h, iz * h };

			r[imorton] = BlockInfo(imorton, idx, origin, h, h / TBlock::sizeX, this->_linaccess(imorton));
		}

		return r;
	}

	public:

	typedef typename TGrid::BlockType BlockType;

	typedef typename TGrid::BlockType TBlock;

	GridMorton(unsigned int nX, unsigned int nY=1, unsigned int nZ=1, const double maxextent = 1):
		TGrid(nX, nY, nZ, maxextent)
	{
		//std::cout << "____________________________________ YOU ARE USING MORTON GRIDS _______________________________" << std::endl;

		_generate_morton_mapping();

		//std::cout << "____________________________________ _generate_morton_mapping _______________________________" << std::endl;

		cached_infos = _getBlocksInfo();
	}

	TBlock& operator()( int ix,  int iy = 0,  int iz = 0) const
	{
		const int N = this->N;
		const int NX = this->NX;
		const int NY = this->NY;
		const int NZ = this->NZ;

		printf("asking for %d %d %d (N: %d %d %d)\n", ix, iy, iz, NX, NY, NZ);
		assert(ix>= -1 && ix <= NX);
		assert(iy>= -1 && iy <= NY);
		assert(iz>= -1 && iz <= NZ);

		const int _ix = (ix + NX) % NX;
		const int _iy = (iy + NY) % NY;
		const int _iz = (iz + NZ) % NZ;

		assert(_ix >= 0 && _ix < NX);
		assert(_iy >= 0 && _iy < NY);
		assert(_iz >= 0 && _iz < NZ);

		const int iflat = _ix + NX * (_iy + NY * _iz);

		assert(iflat >= 0 && iflat < f2m.size());

		const unsigned int idx = f2m[iflat];
		assert(idx >= 0 && idx < N);

		return *(this->_linaccess(idx));
	}

	vector<BlockInfo>& getBlocksInfo()
	{
		return cached_infos;
	}

	const vector<BlockInfo>& getBlocksInfo() const
	{
		return cached_infos;
	}
};
