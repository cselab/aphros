/*
 *  BlockLabMPI.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 10/21/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "GridMPI.h"

template<typename MyBlockLab>
class BlockLabMPI : public MyBlockLab
{
	const SynchronizerMPI * refSynchronizerMPI;
	typedef typename MyBlockLab::BlockType BlockType;
	
protected:
	int mypeindex[3], pesize[3], mybpd[3];
	int gLastX, gLastY, gLastZ;
	
public:
	template< typename TGrid >
	void prepare(GridMPI<TGrid>& grid, const SynchronizerMPI& SynchronizerMPI)
	{
		refSynchronizerMPI = &SynchronizerMPI;
		refSynchronizerMPI->getpedata(mypeindex, pesize, mybpd);
		StencilInfo stencil = refSynchronizerMPI->getstencil();
		assert(stencil.isvalid());
		MyBlockLab::prepare(grid, stencil.sx,  stencil.ex,  stencil.sy,  stencil.ey,  stencil.sz,  stencil.ez, stencil.tensorial);
		gLastX = grid.getBlocksPerDimension(0)-1;
		gLastY = grid.getBlocksPerDimension(1)-1;
		gLastZ = grid.getBlocksPerDimension(2)-1;
	}
	
	void load(const BlockInfo& info, const Real t=0, const bool applybc=true)
	{		
		MyBlockLab::load(info, t, false);
		
		assert(refSynchronizerMPI != NULL);
		
		const int xorigin = mypeindex[0]*mybpd[0];		
		const int yorigin = mypeindex[1]*mybpd[1];
		const int zorigin = mypeindex[2]*mybpd[2];
		
		const bool xskin = (info.index[0] == xorigin || info.index[0] == xorigin + mybpd[0]-1);
		const bool yskin = (info.index[1] == yorigin || info.index[1] == yorigin + mybpd[1]-1);
		const bool zskin = (info.index[2] == zorigin || info.index[2] == zorigin + mybpd[2]-1);
		
		const bool xboundary = info.index[0]==0 || info.index[0]==gLastX;
		const bool yboundary = info.index[1]==0 || info.index[1]==gLastY;
		const bool zboundary = info.index[2]==0 || info.index[2]==gLastZ;		
		
		const bool any_periodic = this->is_xperiodic() || this->is_yperiodic() || this->is_zperiodic();
		const bool any_boundary = xboundary || yboundary || zboundary;
		
		if ((xskin || yskin || zskin))// && (!any_boundary || any_boundary && any_periodic))
		{
			const bool xperiodic = this->is_xperiodic();
			const bool yperiodic = this->is_yperiodic();
			const bool zperiodic = this->is_zperiodic();
			
			const int rsx = (!xperiodic && info.index[0]==0)? 0 : this->m_stencilStart[0];
			const int rex = (!xperiodic && info.index[0]==gLastX) ? BlockType::sizeX : (BlockType::sizeX+this->m_stencilEnd[0]-1);
			const int rsy = (!yperiodic && info.index[1]==0)? 0 : this->m_stencilStart[1];
			const int rey = (!yperiodic && info.index[1]==gLastY) ? BlockType::sizeY : (BlockType::sizeY+this->m_stencilEnd[1]-1);
			const int rsz = (!zperiodic && info.index[2]==0)? 0 : this->m_stencilStart[2];
			const int rez = (!zperiodic && info.index[2]==gLastZ) ? BlockType::sizeZ : (BlockType::sizeZ+this->m_stencilEnd[2]-1);
						
			Real * const dst = (Real *)&this->m_cacheBlock->LinAccess(0);
			
			typedef typename MyBlockLab::ElementType ET;
			
			refSynchronizerMPI->fetch((const Real*)info.ptrBlock, dst, 
									  this->m_stencilStart[0], this->m_stencilStart[1], this->m_stencilStart[2],
									  this->m_cacheBlock->getSize()[0], this->m_cacheBlock->getSize()[1], this->m_cacheBlock->getSize()[2],
									  sizeof(ET)/sizeof(Real),
									  rsx, rex, rsy, rey, rsz, rez);
		}
		
		if (applybc) MyBlockLab::_apply_bc(info, t);
	}
    
    void release()
    {
        MyBlockLab::release();
    }
};
