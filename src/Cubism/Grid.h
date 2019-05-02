/*
 *  Grid.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>

#ifdef _USE_NUMA_
#include <numa.h>
#include <omp.h>
#endif
#include "BlockInfo.h"
#include "MeshMap.h"

//hello git
template <typename B_, template<typename X> class allocator=std::allocator>
class Grid
{
 public:
  using Block = B_;
 
 private:
	Block * m_blocks;
  std::vector<BlockInfo> m_vInfo;

protected:

	const size_t NX, NY, NZ, N;
	const double maxextent;

    const bool m_own_mesh_maps;
    std::vector<MeshMap<Block>*> m_mesh_maps;

	void _dealloc()
	{
		allocator<Block> alloc;

		alloc.deallocate(m_blocks, N);

        if (m_own_mesh_maps)
        {
            for (size_t i = 0; i < m_mesh_maps.size(); ++i)
            {
                delete m_mesh_maps[i];
                m_mesh_maps[i] = NULL;
            }
        }
	}

	void _alloc()
	{
		allocator<Block> alloc;
		m_blocks = alloc.allocate(N);
		assert(m_blocks!=NULL);

		//numa touch
    #if 0
		#pragma omp parallel
		{
#ifdef _USE_NUMA_
			const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
			const int mynode = omp_get_thread_num() / cores_per_node;
			numa_run_on_node(mynode);
#endif
#pragma omp for schedule(static)
			for(int i=0; i<(int)N; ++i)
				m_blocks[i].clear();
		}
    #endif
	}

	Block* _linaccess(const size_t idx) const
	{
		assert(idx >= 0);
		assert(idx < N);

		return m_blocks + idx;
	}

	size_t _encode(const size_t ix, const size_t iy, const size_t iz) const
	{
		assert(ix>=0 && ix<NX);
		assert(iy>=0 && iy<NY);
		assert(iz>=0 && iz<NZ);

		return ix + NX*(iy + NY*iz);
	}

public:

	typedef Block BlockType;

    Grid(const size_t NX, const size_t NY = 1, const size_t NZ = 1, const double maxextent = 1) :
        m_blocks(NULL), NX(NX), NY(NY), NZ(NZ), N(NX*NY*NZ), maxextent(maxextent),
        m_own_mesh_maps(true)
    {
        _alloc();

        // cells per block
        const size_t bs[3] = {Block::bx, Block::by, Block::bz};
        // cells total
        const size_t nn[3] = {NX * bs[0], NY * bs[1], NZ * bs[2]};
        // cell size (h_gridpoint from BlockInfo)
        const double hc = 
            (maxextent / std::max(nn[0], std::max(nn[1], nn[2])));
        // block extent
        const double h = std::max(bs[0], std::max(bs[1], bs[2])) * hc;

        const double extents[3] = {hc*nn[0], hc*nn[1], hc*nn[2]};
        const size_t nBlocks[3] = {NX, NY, NZ};
        for (int i = 0; i < 3; ++i)
        {
            MeshMap<Block>* m = new MeshMap<Block>(
                0.0, extents[i], nBlocks[i], bs[i]);
            UniformDensity uniform;
            m->init(&uniform); // uniform only for this constructor
            m_mesh_maps.push_back(m);
        }

        for(size_t iz=0; iz<NZ; iz++)
            for(size_t iy=0; iy<NY; iy++)
                for(size_t ix=0; ix<NX; ix++)
                {
                    const long long blockID = _encode(ix, iy, iz);
                    const int idx[3] = {(int)ix, (int)iy, (int)iz};
                    const double origin[3] = {ix*h, iy*h, iz*h};

                    m_vInfo.push_back(BlockInfo(blockID, idx, origin, h, hc, _linaccess(blockID)));
                }
    }

    Grid(MeshMap<Block>* const mapX, MeshMap<Block>* const mapY, MeshMap<Block>* const mapZ,
            const int _NX, const int _NY=1, const int _NZ=1) :
        m_blocks(NULL),
        NX(_NX), NY(_NY), NZ(_NZ),
        N(_NX*_NY*_NZ),
        maxextent(-1.0), // not used
        m_own_mesh_maps(false)
    {
        _alloc();

        m_mesh_maps.push_back(mapX);
        m_mesh_maps.push_back(mapY);
        m_mesh_maps.push_back(mapZ);

        for(size_t iz=0; iz<NZ; iz++)
            for(size_t iy=0; iy<NY; iy++)
                for(size_t ix=0; ix<NX; ix++)
                {
                    const long long blockID = _encode(ix, iy, iz);
                    const int idx[3] = {int(ix), int(iy), int(iz)};

                    m_vInfo.push_back(BlockInfo(blockID, idx, mapX, mapY, mapZ, _linaccess(blockID)));
                }
    }

	virtual ~Grid() { _dealloc(); }

	void setup(const size_t nX, const size_t nY, const size_t nZ)
	{
        // TODO: [fabianw@mavt.ethz.ch; Wed May 03 2017 04:30:53 PM (-0700)]
        // dead code (?)
		std::cout << "Setting up the grid with " << nX << "x" << nY << "x" << nZ << " blocks ...";

		_dealloc();

		_alloc();

		std::cout << "done. " << std::endl;
	}

	virtual size_t getBlocksPerDimension(int idim) const
	{
		assert(idim>=0 && idim<3);

		switch (idim)
		{
			case 0: return NX;
			case 1: return NY;
			case 2: return NZ;
			default: abort();
				return 0;
		}
	}

	virtual bool avail(int, int, int) const { return true; }

	virtual Block& operator()(size_t ix, size_t iy=0, size_t iz=0) const
	{
		return *_linaccess( _encode((ix+NX) % NX, (iy+NY) % NY, (iz+NZ) % NZ) );
	}

	virtual std::vector<BlockInfo>& getBlocksInfo()
	{
        return m_vInfo;
	}

	virtual const std::vector<BlockInfo>& getBlocksInfo() const
	{
        return m_vInfo;
	}

    double getH() const
    {
        std::vector<BlockInfo> vInfo = this->getBlocksInfo();
        BlockInfo info = vInfo[0];
        return info.h_gridpoint;
    }

    inline MeshMap<Block>& getMeshMap(const int i)
    {
        assert(i>=0 && i<3);
        return *m_mesh_maps[i];
    }
    inline const MeshMap<Block>& getMeshMap(const int i) const
    {
        assert(i>=0 && i<3);
        return *m_mesh_maps[i];
    }
};

template <typename Block, template<typename X> class allocator>
std::ostream& operator<< (std::ostream& out, const Grid<Block, allocator>& grid)
{
	//save metadata
	out << grid.getBlocksPerDimension(0) << " "
	<< grid.getBlocksPerDimension(1) << " "
	<< grid.getBlocksPerDimension(2) << std::endl;

	return out;
}


template <typename Block, template<typename X> class allocator>
std::ifstream& operator>> (std::ifstream& in, Grid<Block, allocator>& grid)
{
	//read metadata
	size_t nx, ny, nz;
	in >> nx;
    in.ignore(1,' ');
    in >> ny;
    in.ignore(1,' ');
    in >> nz;
	in.ignore(1,'\n');

	grid.setup(nx, ny, nz);

	return in;
}
