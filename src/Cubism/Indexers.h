/*
 *  Indexers.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 10/13/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cassert>
#include <cmath>

class Indexer
{
protected:
	const unsigned int sizeX, sizeY, sizeZ, sizeTotal;
    
public:	
	Indexer(const unsigned int sizeX, const unsigned int sizeY, const unsigned int sizeZ):
    sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), sizeTotal(sizeX*sizeY*sizeZ)
	{
	}
	
	virtual unsigned int encode(unsigned int ix, unsigned int iy, unsigned int iz) const
	{
		const unsigned int retval = ix + sizeX*(iy + iz*sizeY);

		assert(retval < sizeTotal && retval>=0);

		return retval; 
	}
	
	virtual void decode(unsigned int code, unsigned int& ix, unsigned int& iy, unsigned int& iz) const
	{
		ix = code % sizeX;
		iy = (code/sizeX) % sizeY;
		iz = (code/sizeX/sizeY);
	}
};

class IndexerMorton : public Indexer
{
    unsigned int depth;
public:
	IndexerMorton(const unsigned int sizeX, const unsigned int sizeY, const unsigned int sizeZ):
    Indexer(sizeX, sizeY, sizeZ)
	{
        depth = (unsigned int) fmin(10., ceil(log2((double)fmax(sizeX,fmax(sizeY,sizeZ)))));
	}
	
	unsigned int encode(unsigned int ix, unsigned int iy, unsigned int iz) const
	{
        unsigned int idx=0;
        
        for(unsigned int counter=0;counter<depth;++counter)
        {
            const unsigned int bitmask = 1 << counter;
            const unsigned int idx0 = ix&bitmask;
            const unsigned int idx1 = iy&bitmask;
            const unsigned int idx2 = iz&bitmask;
            
            idx |= ((idx0<<2*counter) | (idx1<<(2*counter+1)) | (idx2<<(2*counter+2)));
        }
        
		return idx; 
	}
	
	void decode(unsigned int code, unsigned int& ix, unsigned int& iy, unsigned int& iz) const
	{
        ix=iy=iz=0;
        
        for(unsigned int counter=0;counter<depth;++counter)
        {
            const unsigned int bitmask_x = 1 << (counter*3+0);
            const unsigned int bitmask_y = 1 << (counter*3+1);
            const unsigned int bitmask_z = 1 << (counter*3+2);
            ix |= (code&bitmask_x)>>2*counter;
            iy |= (code&bitmask_y)>>(2*counter+1);
            iz |= (code&bitmask_z)>>(2*counter+2);
        }
	}
};
