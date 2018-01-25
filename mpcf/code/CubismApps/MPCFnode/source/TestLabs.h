/*
 *  TestLabs.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger on 10/20/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TESTLABS_H_DQM59CVK
#define TESTLABS_H_DQM59CVK

#include <cassert>
#include <vector>
#include <string>
#include <cmath>
#include "Types.h"
#include "BlockLab.h"
#include "BoundaryConditions.h"


///////////////////////////////////////////////////////////////////////////////
// Helper
///////////////////////////////////////////////////////////////////////////////
template <typename TBlock>
static void _assign_grid_spacing(const BlockInfo& info, Real h[6]
#ifndef NDEBUG
        , const bool check[6]
#endif /* NDEBUG */
        )
{
    const int N[3] = {TBlock::sizeX-1, TBlock::sizeY-1, TBlock::sizeZ-1};
    for (int j = 0; j < 3; ++j)
    {
        if (info.bUniform[j])
        {
            h[0+2*j] = info.uniform_grid_spacing[j];
            h[1+2*j] = info.uniform_grid_spacing[j];
        }
        else
        {
            const double* const delta = info.ptr_grid_spacing[j];
            h[0+2*j] = delta[0];
            h[1+2*j] = delta[N[j]];
#ifndef NDEBUG
            if (check[0+2*j])
            {
                for (int i = 0; i < 6; ++i)
                    assert(h[0+2*j]==delta[i]);
            }
            if (check[1+2*j])
            {
                for (int i = 0; i < 6; ++i)
                    assert(h[1+2*j]==delta[N[j]-i]);
            }
#endif /* NDEBUG */
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// Cloud cases
///////////////////////////////////////////////////////////////////////////////

template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloud1DCharNonReflectAcousticForcing_5eq: public BlockLab<BlockType,Alloc>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:

    virtual inline std::string name() const { return "BlockLabCloud_1DCharNonReflect_AcousticForcing"; }
    bool is_xperiodic() {return true;}
    bool is_yperiodic() {return true;}
    bool is_zperiodic() {return true;}

    BlockLabCloud1DCharNonReflectAcousticForcing_5eq()
      : BlockLab<BlockType,Alloc>(){}
};


#endif /* TESTLABS_H_DQM59CVK */
