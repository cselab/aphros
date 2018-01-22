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

// CAUTION
// for pressure relaxation: Poisson equation is solved and pressure is assumed
// to be stored in the energy component
template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudLaplace_5eq: public BlockLab<BlockType,Alloc>
{
	typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloudLaplace_5eq"; }
	bool is_xperiodic() {return false;}
	bool is_yperiodic() {return false;}
	bool is_zperiodic() {return false;}

	BlockLabCloudLaplace_5eq(): BlockLab<BlockType,Alloc>(){}

	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        /* TODO: (Thu 07 May 2015 05:28:19 PM CEST) dangerous to hard-code this
         * boundary data! */
        const Real G = 1.0/(Simulation_Environment::GAMMA1-1.0);
        const Real F = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
        ElementTypeBlock b;
        b.clear();
        b.alpha1rho1 = CloudData::rho1;
        b.alpha2rho2 = 0.0;
        b.ru = b.rv = b.rw = 0;
        b.alpha2 = 0;
        b.energy = CloudData::p1;

        if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(b);
        if (info.index[0]==this->NX-1)  bc.template applyBC_dirichlet<0,1>(b);
        if (info.index[1]==0)	        bc.template applyBC_dirichlet<1,0>(b);
        if (info.index[1]==this->NY-1)	bc.template applyBC_dirichlet<1,1>(b);
        if (info.index[2]==0)		bc.template applyBC_dirichlet<2,0>(b);
        if (info.index[2]==this->NZ-1)	bc.template applyBC_dirichlet<2,1>(b);
    }
};

template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudDBC_5eq: public BlockLab<BlockType,Alloc>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloudDBC_5eq"; }
        bool is_xperiodic() {return false;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return false;}

        BlockLabCloudDBC_5eq(): BlockLab<BlockType,Alloc>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        /* TODO: (Thu 07 May 2015 05:28:19 PM CEST) dangerous to hard-code this
 *          * boundary data! */
        const Real G = 1.0/(Simulation_Environment::GAMMA1-1.0);
        const Real F = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
        ElementTypeBlock b;
        b.clear();
        b.alpha1rho1 = CloudData::rho1;
        b.alpha2rho2 = 0.0;
        b.ru = b.rv = b.rw = 0;
        b.alpha2 = 0;
        b.energy = CloudData::p1*G + F*G;

        if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(b);
        if (info.index[0]==this->NX-1)  bc.template applyBC_dirichlet<0,1>(b);
        if (info.index[1]==0)           bc.template applyBC_dirichlet<1,0>(b);
        if (info.index[1]==this->NY-1)  bc.template applyBC_dirichlet<1,1>(b);
        if (info.index[2]==0)           bc.template applyBC_dirichlet<2,0>(b);
        if (info.index[2]==this->NZ-1)  bc.template applyBC_dirichlet<2,1>(b);
    }
};

template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudFarField_5eq: public BlockLab<BlockType,Alloc>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloudFarField_5eq"; }
        bool is_xperiodic() {return false;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return false;}

        BlockLabCloudFarField_5eq(): BlockLab<BlockType,Alloc>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        ElementTypeBlock b;
        // note: primitives are assumed for the reference element (far-field state)
        b.clear();
        b.alpha1rho1 = CloudData::rho1;
        b.alpha2rho2 = 0.0;
        b.ru = b.rv = b.rw = 0.0;
        b.alpha2 = 0.0;
        b.energy = CloudData::p1;

#ifndef _FFEXTRA_
        if (info.index[0]==0)           bc.template applyBC_farfield<0,0>(b);
        if (info.index[0]==this->NX-1)  bc.template applyBC_farfield<0,1>(b);
        if (info.index[1]==0)           bc.template applyBC_farfield<1,0>(b);
        if (info.index[1]==this->NY-1)  bc.template applyBC_farfield<1,1>(b);
        if (info.index[2]==0)           bc.template applyBC_farfield<2,0>(b);
        if (info.index[2]==this->NZ-1)  bc.template applyBC_farfield<2,1>(b);
#else
        if (info.index[0]==0)           bc.template applyBC_farfield_extra<0,0>(b);
        if (info.index[0]==this->NX-1)  bc.template applyBC_farfield_extra<0,1>(b);
        if (info.index[1]==0)           bc.template applyBC_farfield_extra<1,0>(b);
        if (info.index[1]==this->NY-1)  bc.template applyBC_farfield_extra<1,1>(b);
        if (info.index[2]==0)           bc.template applyBC_farfield_extra<2,0>(b);
        if (info.index[2]==this->NZ-1)  bc.template applyBC_farfield_extra<2,1>(b);
#endif
    }
};


template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudDirichletAcousticForcing_5eq: public BlockLab<BlockType,Alloc>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:

    virtual inline std::string name() const { return "BlockLabCloud_Dirichlet_AcousticForcing"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    BlockLabCloudDirichletAcousticForcing_5eq(): BlockLab<BlockType,Alloc>(){}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        ElementTypeBlock myElement = CloudDataAcoustic::boundaryElement0;
        const double pressure = CloudData::p1 + CloudDataAcoustic::p_amplitude*std::sin(2.0*M_PI*CloudDataAcoustic::frequency * t);

        const double G = 1.0/(Simulation_Environment::GAMMA1 - 1.0);
        const double F = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
        myElement.energy = G*pressure + F*G; // zero velocity boundary!

        if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(myElement);
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
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

        //if (info.index[0]==0)
        //{
        //    BGC::BoundaryBlock& myBblock = BGC::bgblock_list_dir0_side0[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)];
        //    for (int i = 0; i < 3*_BLOCKSIZE_*_BLOCKSIZE_; ++i)
        //        myBblock.fe[i].pressure = pressure;
        //
        //    bc.template applyBC_1dchar_ode<0,0>(myBblock);
        //}

        //if (info.index[0]==0)           bc.template applyBC_absorbing<0,0>();
        //if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();

        if (info.index[0]==0)           bc.template applyBC_1dchar_ode<0,0>(BGC::bgblock_list_dir0_side0[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)]);
        if (info.index[0]==this->NX-1)  bc.template applyBC_1dchar_ode<0,1>(BGC::bgblock_list_dir0_side1[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)]);

        #ifdef _NR_CHAR_Y_
        if (info.index[1]==0)           bc.template applyBC_1dchar_ode<1,0>(BGC::bgblock_list_dir1_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
        if (info.index[1]==this->NY-1)  bc.template applyBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
        #else
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        #endif

        //if (info.index[1]==0)
        //{
        //  bc.template applyBC_noslip_simple<1,0>();
        //  bc.template applyBC_noslip_simple_tensorial_pbc<1,0>();
        //}
        //if (info.index[1]==this->NY-1)
        //{
        //  bc.template applyBC_noslip_simple<1,1>();
        //  bc.template applyBC_noslip_simple_tensorial_pbc<1,1>();
        //}
        
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

        //if (info.index[0]==0)           bc.template updateBC_1dchar_ode<0,0>(BGC::bgblock_list_dir0_side0[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)],
        //        h[0], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
        //if (info.index[0]==this->NX-1)  bc.template updateBC_1dchar_ode<0,1>(BGC::bgblock_list_dir0_side1[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)],
        //        h[1], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
        // if (info.index[1]==0)           bc.template updateBC_1dchar_ode<1,0>(BGC::bgblock_list_dir1_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)],
        //         h[2], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
        // if (info.index[1]==this->NY-1)  bc.template updateBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)],
        //         h[3], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
        // if (info.index[2]==0)           bc.template updateBC_1dchar_ode<2,0>(BGC::bgblock_list_dir2_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)],
        //         h[4], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
        // if (info.index[2]==this->NZ-1)  bc.template updateBC_1dchar_ode<2,1>(BGC::bgblock_list_dir2_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)],
        //         h[5], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
    }
};


template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloud1DCharNonReflect_5eq: public BlockLab<BlockType,Alloc>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloud1DCharNonReflect_5eq"; }
        bool is_xperiodic() {return false;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return false;}

        BlockLabCloud1DCharNonReflect_5eq(): BlockLab<BlockType,Alloc>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
          BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

          if (info.index[0]==0)           bc.template applyBC_1dchar_ode<0,0>(BGC::bgblock_list_dir0_side0[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)]);
          if (info.index[0]==this->NX-1)  bc.template applyBC_1dchar_ode<0,1>(BGC::bgblock_list_dir0_side1[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)]);
          if (info.index[1]==0)           bc.template applyBC_1dchar_ode<1,0>(BGC::bgblock_list_dir1_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
          if (info.index[1]==this->NY-1)  bc.template applyBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
          if (info.index[2]==0)           bc.template applyBC_1dchar_ode<2,0>(BGC::bgblock_list_dir2_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)]);
          if (info.index[2]==this->NZ-1)  bc.template applyBC_1dchar_ode<2,1>(BGC::bgblock_list_dir2_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)]);
        }

        void apply_bc_update(const BlockInfo& info, const Real dt=0, const Real a=0, const Real b=0)
        {
          BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

          Real h[6];
#ifndef NDEBUG
        const bool check[6] = {true};
#endif /* NDEBUG */
        _assign_grid_spacing<BlockType>(info, h
#ifndef NDEBUG
                , check
#endif /* NDEBUG */
                );

          if (info.index[0]==0)           bc.template updateBC_1dchar_ode<0,0>(BGC::bgblock_list_dir0_side0[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)],
                  h[0], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[0]==this->NX-1)  bc.template updateBC_1dchar_ode<0,1>(BGC::bgblock_list_dir0_side1[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)],
                  h[1], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[1]==0)           bc.template updateBC_1dchar_ode<1,0>(BGC::bgblock_list_dir1_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)],
                  h[2], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[1]==this->NY-1)  bc.template updateBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)],
                  h[3], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[2]==0)           bc.template updateBC_1dchar_ode<2,0>(BGC::bgblock_list_dir2_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)],
                  h[4], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[2]==this->NZ-1)  bc.template updateBC_1dchar_ode<2,1>(BGC::bgblock_list_dir2_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)],
                  h[5], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
        }
};

template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudSym1DCharNonReflect_5eq: public BlockLab<BlockType,Alloc>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloudSym1DCharNonReflect_5eq"; }
        bool is_xperiodic() {return false;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return false;}

        BlockLabCloudSym1DCharNonReflect_5eq(): BlockLab<BlockType,Alloc>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
          BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

          if (info.index[0]==0)           bc.template applyBC_symmetry<0,0>();
          if (info.index[0]==this->NX-1)  bc.template applyBC_1dchar_ode<0,1>(BGC::bgblock_list_dir0_side1[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)]);
          if (info.index[1]==0)           bc.template applyBC_symmetry<1,0>();
          if (info.index[1]==this->NY-1)  bc.template applyBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
          if (info.index[2]==0)           bc.template applyBC_symmetry<2,0>();
          if (info.index[2]==this->NZ-1)  bc.template applyBC_1dchar_ode<2,1>(BGC::bgblock_list_dir2_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)]);
        }

        void apply_bc_update(const BlockInfo& info, const Real dt=0, const Real a=0, const Real b=0)
        {
          BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

          Real h[6];
#ifndef NDEBUG
        const bool check[6] = {false, true, false, true, false, true};
#endif /* NDEBUG */
        _assign_grid_spacing<BlockType>(info, h
#ifndef NDEBUG
                , check
#endif /* NDEBUG */
                );

          if (info.index[0]==this->NX-1)  bc.template updateBC_1dchar_ode<0,1>(BGC::bgblock_list_dir0_side1[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)],
                  h[1], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[1]==this->NY-1)  bc.template updateBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)],
                  h[3], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[2]==this->NZ-1)  bc.template updateBC_1dchar_ode<2,1>(BGC::bgblock_list_dir2_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)],
                  h[5], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
        }
};


template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudWall1DCharNonReflect_5eq: public BlockLab<BlockType,Alloc>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloudWall1DCharNonReflect_5eq"; }
        bool is_xperiodic() {return false;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return false;}

        BlockLabCloudWall1DCharNonReflect_5eq(): BlockLab<BlockType,Alloc>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
          BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

          if (info.index[0]==0)           bc.template applyBC_1dchar_ode<0,0>(BGC::bgblock_list_dir0_side0[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)]);
          if (info.index[0]==this->NX-1)  bc.template applyBC_symmetry<0,1>();
          if (info.index[1]==0)           bc.template applyBC_1dchar_ode<1,0>(BGC::bgblock_list_dir1_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
          if (info.index[1]==this->NY-1)  bc.template applyBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)]);
          if (info.index[2]==0)           bc.template applyBC_1dchar_ode<2,0>(BGC::bgblock_list_dir2_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)]);
          if (info.index[2]==this->NZ-1)  bc.template applyBC_1dchar_ode<2,1>(BGC::bgblock_list_dir2_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)]);
        }

        void apply_bc_update(const BlockInfo& info, const Real dt=0, const Real a=0, const Real b=0)
        {
          BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

          Real h[6];
#ifndef NDEBUG
        const bool check[6] = {true, false, true, true, true, true};
#endif /* NDEBUG */
        _assign_grid_spacing<BlockType>(info, h
#ifndef NDEBUG
                , check
#endif /* NDEBUG */
                );

          if (info.index[0]==0)           bc.template updateBC_1dchar_ode<0,0>(BGC::bgblock_list_dir0_side0[info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ)],
                  h[0], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[1]==0)           bc.template updateBC_1dchar_ode<1,0>(BGC::bgblock_list_dir1_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)],
                  h[2], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[1]==this->NY-1)  bc.template updateBC_1dchar_ode<1,1>(BGC::bgblock_list_dir1_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ)],
                  h[3], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[2]==0)           bc.template updateBC_1dchar_ode<2,0>(BGC::bgblock_list_dir2_side0[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)],
                  h[4], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
          if (info.index[2]==this->NZ-1)  bc.template updateBC_1dchar_ode<2,1>(BGC::bgblock_list_dir2_side1[info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY)],
                  h[5], dt, a, b, BGC::pamb, BGC::L, BGC::lambda);
        }
};


template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudAbsorb_5eq: public BlockLab<BlockType,Alloc>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloudAbsorb_5eq"; }
        bool is_xperiodic() {return false;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return false;}

        BlockLabCloudAbsorb_5eq(): BlockLab<BlockType,Alloc>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)
        {
            bc.template applyBC_absorbing<0,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<0,0>();
        }
        if (info.index[0]==this->NX-1)
        {
            bc.template applyBC_absorbing<0,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<0,1>();
        }
        if (info.index[1]==0)
        {
            bc.template applyBC_absorbing<1,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<1,0>();
        }
        if (info.index[1]==this->NY-1)
        {
            bc.template applyBC_absorbing<1,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<1,1>();
        }
        if (info.index[2]==0)
        {
            bc.template applyBC_absorbing<2,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<2,0>();
        }
        if (info.index[2]==this->NZ-1)
        {
            bc.template applyBC_absorbing<2,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<2,1>();
        }
    }
};


template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudSymAbsorb_5eq: public BlockLab<BlockType,Alloc>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloudSymAbsorb_5eq"; }
        bool is_xperiodic() {return false;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return false;}

        BlockLabCloudSymAbsorb_5eq(): BlockLab<BlockType,Alloc>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
/*
        if (info.index[0]==0)
            bc.template applyBC_symmetry<0,0>();
        if (info.index[0]==this->NX-1)
            bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)
            bc.template applyBC_symmetry<1,0>();
        if (info.index[1]==this->NY-1)
            bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)
            bc.template applyBC_symmetry<2,0>();
        if (info.index[2]==this->NZ-1)
            bc.template applyBC_absorbing<2,1>();
*/

        if (info.index[0]==0)
        {
            bc.template applyBC_symmetry<0,0>();
            bc.template applyBC_symmetry_tensorials_edges_corner<0,0>();
        }
        if (info.index[0]==this->NX-1)
        {
            bc.template applyBC_absorbing<0,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<0,1>();
        }
        if (info.index[1]==0)
        {
            bc.template applyBC_symmetry<1,0>();
            bc.template applyBC_symmetry_tensorials_edges_corner<1,0>();
        }
        if (info.index[1]==this->NY-1)
        {
            bc.template applyBC_absorbing<1,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<1,1>();
        }
        if (info.index[2]==0)
        {
            bc.template applyBC_symmetry<2,0>();
            bc.template applyBC_symmetry_tensorials_edges_corner<2,0>();
        }
        if (info.index[2]==this->NZ-1)
        {
            bc.template applyBC_absorbing<2,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<2,1>();
        }


// basically the same: tensorial kernels not considered
/*
        if (info.index[0]==0)
            bc.template applyBC_reflecting<0,0>();
        if (info.index[0]==this->NX-1)
            bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)
            bc.template applyBC_reflecting<1,0>();
        if (info.index[1]==this->NY-1)
            bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)
            bc.template applyBC_reflecting<2,0>();
        if (info.index[2]==this->NZ-1)
            bc.template applyBC_absorbing<2,1>();
*/
    }
};


template<typename BlockType, template<typename X> class Alloc=std::allocator>
class BlockLabCloudWallAbsorb_5eq: public BlockLab<BlockType,Alloc>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:

        virtual inline std::string name() const { return "BlockLabCloudWallAbsorb_5eq"; }
        bool is_xperiodic() {return false;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return false;}

        BlockLabCloudWallAbsorb_5eq(): BlockLab<BlockType,Alloc>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
        BoundaryCondition<BlockType,ElementTypeBlock,Alloc> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)
        {
            bc.template applyBC_absorbing<0,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<0,0>();
        }
        if (info.index[0]==this->NX-1)
        {
            bc.template applyBC_symmetry<0,1>();
            bc.template applyBC_symmetry_tensorials_edges_corner<0,1>();
        }
        if (info.index[1]==0)
        {
            bc.template applyBC_absorbing<1,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<1,0>();
        }
        if (info.index[1]==this->NY-1)
        {
            bc.template applyBC_absorbing<1,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<1,1>();
        }
        if (info.index[2]==0)
        {
            bc.template applyBC_absorbing<2,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<2,0>();
        }
        if (info.index[2]==this->NZ-1)
        {
            bc.template applyBC_absorbing<2,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<2,1>();
        }
       }
};


///////////////////////////////////////////////////////////////////////////////
// SIC cases
///////////////////////////////////////////////////////////////////////////////
struct SICCloudBCData
{
    static FluidElement boundaryElement0;
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabSICCloud: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    BlockLabSICCloud(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "BlockLabSICCloud"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)           bc.template applyBC_absorbing<0,0>();
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabSICCloud_AbsorbingTensorial: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    BlockLabSICCloud_AbsorbingTensorial(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "BlockLabSICCloud_AbsorbingTensorial"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)
        {
            bc.template applyBC_absorbing<0,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<0,0>();
        }
        if (info.index[0]==this->NX-1)
        {
            bc.template applyBC_absorbing<0,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<0,1>();
        }
        if (info.index[1]==0)
        {
            bc.template applyBC_absorbing<1,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<1,0>();
        }
        if (info.index[1]==this->NY-1)
        {
            bc.template applyBC_absorbing<1,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<1,1>();
        }
        if (info.index[2]==0)
        {
            bc.template applyBC_absorbing<2,0>();
            bc.template applyBC_absorbing_tensorials_edges_corner<2,0>();
        }
        if (info.index[2]==this->NZ-1)
        {
            bc.template applyBC_absorbing<2,1>();
            bc.template applyBC_absorbing_tensorials_edges_corner<2,1>();
        }
    }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabSICCloud_Dirichlet: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    BlockLabSICCloud_Dirichlet(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "BlockLabSICCloud_Dirichlet"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(SICCloudBCData::boundaryElement0);
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabSICCloud_LeftReflect: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    BlockLabSICCloud_LeftReflect(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "BlockLabSICCloud_LeftReflect"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)           bc.template applyBC_reflecting<0,0>();
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        // if (info.index[0]==this->NX-1)  bc.template applyBC_dirichlet<0,1>(SICCloudBCData::boundaryElement0);
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
};


///////////////////////////////////////////////////////////////////////////////
// Channel flow cases
///////////////////////////////////////////////////////////////////////////////
template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabChannelSimple: public BlockLab<BlockType,allocator>
{
        typedef typename BlockType::ElementType ElementTypeBlock;

public:
        virtual std::string name() const { return "BlockLabChannel"; }
        bool is_xperiodic() {return true;}
        bool is_yperiodic() {return false;}
        bool is_zperiodic() {return true;}

        BlockLabChannelSimple(): BlockLab<BlockType,allocator>(){}

        void _apply_bc(const BlockInfo& info, const Real t=0)
        {
          BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

          if (info.index[1]==0)
          {
            bc.template applyBC_noslip_simple<1,0>();
            bc.template applyBC_noslip_simple_tensorial_pbc<1,0>();
          }
          if (info.index[1]==this->NY-1)
          {
            bc.template applyBC_noslip_simple<1,1>();
            bc.template applyBC_noslip_simple_tensorial_pbc<1,1>();
          }
        }
};

///////////////////////////////////////////////////////////////////////////////
// Advection cases
///////////////////////////////////////////////////////////////////////////////
struct AdvectionBCData
{
    static std::vector<FluidElement> boundaryElementVector;
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabAdvection_DirichletUniform: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    BlockLabAdvection_DirichletUniform(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "BlockLabAdvection_DirichletUniform"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(AdvectionBCData::boundaryElementVector[0]);
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabAdvection_DirichletNonUniform: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    BlockLabAdvection_DirichletNonUniform(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "BlockLabAdvection_DirichletNonUniform"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)           bc.template applyBC_dirichletNonUniform<0,0>(info, AdvectionBCData::boundaryElementVector);
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
};


#endif /* TESTLABS_H_DQM59CVK */
