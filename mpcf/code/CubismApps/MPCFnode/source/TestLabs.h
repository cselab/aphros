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
#include "NonUniform.h"
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
namespace CloudData
{
    extern double rho1, rho2;
    extern double p1, p2;
}

struct CloudDataAcoustic
{
    static FluidElement boundaryElement0;
    static double p_amplitude;
    static double p_amplitudea;
    static double frequency;
    static double t0;
    static double phase0;
    static double sigma;
};


template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloud1DCharNonReflectAcousticForcing_5eq: public BlockLab<BlockType,Alloc>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:

    virtual inline std::string name() const { return "BlockLabCloud_1DCharNonReflect_AcousticForcing"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    #ifdef _2D_
    bool is_zperiodic() {return true;}
    #else
    bool is_zperiodic() {return false;}
    #endif

    BlockLabCloud1DCharNonReflectAcousticForcing_5eq(): BlockLab<BlockType,Alloc>(){}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        // Gaussian modulated sinusoid:
        // exp(-((t-t0) * freq / sigma) ** 2) * cos(2. * pi * (freq * (t-t0) - phase0))
        typedef CloudDataAcoustic Cda;
        const double pressure = CloudData::p1 + Cda::p_amplitude
          * std::exp(-std::pow((t - Cda::t0) * Cda::frequency / Cda::sigma, 2))
          * std::cos(2.0 * M_PI * (Cda::frequency * (t - Cda::t0) - Cda::phase0));

        if (info.index[0]==0)           bc.template applyBC_1dchar_ode<0,0>(BGC::bgblock_list_dir0_side0[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)]);
        if (info.index[0]==this->NX-1)  bc.template applyBC_1dchar_ode<0,1>(BGC::bgblock_list_dir0_side1[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)]);

        #ifdef _NR_CHAR_Y_
        if (info.index[1]==0)           bc.template applyBC_1dchar_ode<1,0>(BGC::bgblock_list_dir1_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
        if (info.index[1]==this->NY-1)  bc.template applyBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
        #else
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        #endif
        
        #ifndef _2D_
        #ifdef _NR_CHAR_Z_
        if (info.index[2]==0)           bc.template applyBC_1dchar_ode<2,0>(BGC::bgblock_list_dir2_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)]);
        if (info.index[2]==this->NZ-1)  bc.template applyBC_1dchar_ode<2,1>(BGC::bgblock_list_dir2_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)]);
        #else
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
        #endif
        #endif
    }

    void apply_bc_update(const BlockInfo& info, const Real dt=0, const Real a=0, const Real b=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        Real h[6];
#ifndef NDEBUG
        const bool check[6] = {true, true, false, false, false, false};
#endif /* NDEBUG */
        _assign_grid_spacing<BlockType>(info, h
#ifndef NDEBUG
                , check
#endif /* NDEBUG */
                );

    }
};


#endif /* TESTLABS_H_DQM59CVK */
