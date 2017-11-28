// File       : BoundaryConditions.h
// Date       : Fri 01 Apr 2016 05:58:23 PM CEST
// Author     : Fabian Wermelinger
// Description: Some common boundary conditions
// Copyright 2016 ETH Zurich. All Rights Reserved.
#ifndef BOUNDARYCONDITIONS_H_OQJPQNBH
#define BOUNDARYCONDITIONS_H_OQJPQNBH

#include "Types.h"
#include "Matrix3D.h"

template<typename TBlock, typename TElement, template<typename X> class allocator=std::allocator>
class BoundaryCondition
{
protected:

    int s[3], e[3];
    int stencilStart[3], stencilEnd[3];
    Matrix3D<TElement, true, allocator> * cacheBlock;

    template<int dir, int side>
    void _setup()
    {
        s[0] =  dir==0? (side==0? stencilStart[0]: TBlock::sizeX) : 0;
        s[1] =  dir==1? (side==0? stencilStart[1]: TBlock::sizeY) : 0;
        s[2] =  dir==2? (side==0? stencilStart[2]: TBlock::sizeZ) : 0;

        e[0] =  dir==0? (side==0? 0: TBlock::sizeX + stencilEnd[0]-1) : TBlock::sizeX;
        e[1] =  dir==1? (side==0? 0: TBlock::sizeY + stencilEnd[1]-1) : TBlock::sizeY;
        e[2] =  dir==2? (side==0? 0: TBlock::sizeZ + stencilEnd[2]-1) : TBlock::sizeZ;
    }

public:

    BoundaryCondition(const int ss[3], const int se[3], Matrix3D<TElement, true, allocator> * cacheBlock):
    cacheBlock(cacheBlock)
    {
        s[0]=s[1]=s[2]=0;
        e[0]=e[1]=e[2]=0;

        stencilStart[0] = ss[0];
        stencilStart[1] = ss[1];
        stencilStart[2] = ss[2];

        stencilEnd[0] = se[0];
        stencilEnd[1] = se[1];
        stencilEnd[2] = se[2];
    }

    TElement& operator()(int ix, int iy, int iz)
    {
        // cacheBlock points to the cacheBlock of the Lab
        return cacheBlock->Access(ix-stencilStart[0],iy-stencilStart[1],iz-stencilStart[2]);
    }

    // 0-th order extrapolation absorbing boundary, used in Types.h for the
    // specific Lab-Type of the application (MyLabAbsorbing)
    template<int dir, int side>
    void applyBC_absorbing()
    {
        _setup<dir,side>();

        for(int iz=s[2]; iz<e[2]; iz++)
            for(int iy=s[1]; iy<e[1]; iy++)
                for(int ix=s[0]; ix<e[0]; ix++)
                {
                    (*this)(ix,iy,iz) = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
                                                dir==1? (side==0? 0:TBlock::sizeY-1):iy,
                                                dir==2? (side==0? 0:TBlock::sizeZ-1):iz);
                }
    }
};

#endif /* BOUNDARYCONDITIONS_H_OQJPQNBH */
