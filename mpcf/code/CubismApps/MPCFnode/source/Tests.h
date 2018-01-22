/*
 *  Tests.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/20/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "TestLabs.h"

/* TODO: (fabianw; Wed 03 Jun 2015 07:28:43 PM CEST) For now, the specific Lab
 * you want to use for your test case MUST be specified here.  We use runtime
 * selection of the test case in the main.cpp file, but Lab (which depends on
 * the chosen test case at runtime) is compile time.  This requires a better
 * solution. */

/* #include "Test_ShockTube.h" */
/* typedef BlockLab<Block_t, std::allocator> Lab; */

//#include "Test_SIC.h"
//maybe replace it with std::allocator
//typedef BlockLabCollapse<Block_t, tbb::scalable_allocator> Lab;
//typedef BlockLabCollapse<Block_t, std::allocator> Lab;

//typedef BlockLabBubbleYZSymmetric<Block_t, tbb::scalable_allocator> Lab;

//#include "Test_DoubleMachReflection.h"
//typedef BlockLabDMR<Block_t, tbb::scalable_allocator> Lab;

//#include "Test_CVT.h"
//typedef BlockLab_CVT<Block_t, tbb::scalable_allocator> Lab;

//#include "Test_SICCloud.h"
//typedef BlockLabSICCloud<Block_t, std::allocator> Lab;

#if defined(_BCLABSODSTD_)
#include "Test_SOD.h"
// standard
typedef BlockLabSOD<Block_t, std::allocator> Lab;

#elif defined(_BCLABSODSHARP_)
#include "Test_SOD.h"
// tensorial, if sharpening is included
typedef BlockLabSODSharp<Block_t, std::allocator> Lab;

#elif defined(_BCLABSODXPER_)
#include "Test_SOD.h"
// tensorial, if sharpening is included
typedef BlockLabSOD_xperiodic<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDLAPLACE_)
// boundary conditions if pressure relaxation is included
/* typedef BlockLabCloudLaplace<Block_t, std::allocator> Lab; */
typedef BlockLabCloudLaplace_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDFARFIELD_)
// far-field boundary condition
typedef BlockLabCloudFarField_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDDIRICHLET_ACOUSTIC_)
// 1d characteristic-based non-refelcting boundary condition with acoustic
// forcing pressure
typedef BlockLabCloudDirichletAcousticForcing_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUD1DCHARNREF_ACOUSTIC_)
// 1d characteristic-based non-refelcting boundary condition with acoustic
// forcing pressure
#define _CHARACTERISTIC_1D_BOUNDARY_
typedef BlockLabCloud1DCharNonReflectAcousticForcing_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUD1DCHARNREF_)
// 1d characteristic-based non-refelcting boundary condition
#define _CHARACTERISTIC_1D_BOUNDARY_
typedef BlockLabCloud1DCharNonReflect_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDSYM1DCHARNREF_)
// 1d characteristic-based non-refelcting boundary condition at 3 sides and symmetry boundary condition for the remaining ones
#define _CHARACTERISTIC_1D_BOUNDARY_
typedef BlockLabCloudSym1DCharNonReflect_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDDBC_)
// Dirichlet boundary condition
typedef BlockLabCloudDBC_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDABSORB_)
// absorbing boundary condition (tensorial)
typedef BlockLabCloudAbsorb_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDSYMABSORB_)
// symmetry boundary condition at 3 sides combined with absorbing boundary conditions for the remaining ones
typedef BlockLabCloudSymAbsorb_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDWALLABSORB_)
// symmetry/reflecting wall boundary condition at 1 side combined with absorbing boundary conditions for the remaining ones
typedef BlockLabCloudWallAbsorb_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABCLOUDWALL1DCHARNREF_)
// symmetry/reflecting wall boundary condition at 1 side combined with 1d characteristic-based non-refelcting boundary conditions for the remaining ones
#define _CHARACTERISTIC_1D_BOUNDARY_
typedef BlockLabCloudWall1DCharNonReflect_5eq<Block_t, std::allocator> Lab;

#elif defined(_BCLABSB_)
#include "Test_ShockBubble.h"
typedef BlockLabBubble<Block_t, std::allocator> Lab;

#elif defined(_BCLABSICCLOUD_)
typedef BlockLabSICCloud<Block_t, std::allocator> Lab;

#elif defined(_BCLABSICCLOUD_ABSORB_TENSORIAL_)
typedef BlockLabSICCloud_AbsorbingTensorial<Block_t, std::allocator> Lab;

#elif defined(_BCLABSICCLOUD_DIRICHLET_)
typedef BlockLabSICCloud_Dirichlet<Block_t, std::allocator> Lab;

#elif defined(_BCLABSICCLOUD_LEFTREFLECT_)
typedef BlockLabSICCloud_LeftReflect<Block_t, std::allocator> Lab;

// for oscillating shape test (testing surface tension)
#elif defined(_BCLABOSC3D_)
#include "Test_OscillatingShape.h"
typedef BlockLabOscillatingShape<Block_t, std::allocator> Lab;

#elif defined(_BCLABOSC2D_)
#include "Test_OscillatingShape.h"
typedef BlockLabOscillatingShape2D<Block_t, std::allocator> Lab;


#elif defined(_BCLABCHASIMPLE_)
typedef BlockLabChannelSimple<Block_t, std::allocator> Lab;


#elif defined(_BCLABADVECTION_UNIFORM_)
typedef BlockLabAdvection_DirichletUniform<Block_t, std::allocator> Lab;

#elif defined(_BCLABADVECTION_NONUNIFORM_)
typedef BlockLabAdvection_DirichletNonUniform<Block_t, std::allocator> Lab;


#elif defined(_BCLABSHEAR_)
#include "Test_ShearFlow.h"
typedef BlockLabShear<Block_t, std::allocator> Lab;

#elif defined(_TESTLAB_)
template<typename BlockType, template<typename X> class allocator=std::allocator>
class TestLab: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    TestLab(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "TestLab"; }
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
typedef TestLab<Block_t, std::allocator> Lab;

#else
// generic
typedef BlockLab<Block_t, std::allocator> Lab;
#endif

/* // for 2DSBI */
/* #include "Test_2DSBI.h" */
/* /1* typedef BlockLabCloudLaplace<Block_t, std::allocator> Lab; *1/ */
/* typedef BlockLab2DSBI_5eq<Block_t, std::allocator> Lab; */

///////////////////////////////////////////////////////////////////////////////
// TESTING STUFF
///////////////////////////////////////////////////////////////////////////////
// for HIT cloud test special
/* #include "Test_Cloud.h" */
/* template<typename BlockType, template<typename X> class allocator=std::allocator> */
/* class BlockLabHITCloud_special: public BlockLab<BlockType,allocator> */
/* { */
/* 	typedef typename BlockType::ElementType ElementTypeBlock; */

/* public: */

/*     virtual inline std::string name() const { return "BlockLabHITCloud_special"; } */
/* 	bool is_xperiodic() {return false;} */
/* 	bool is_yperiodic() {return false;} */
/* 	bool is_zperiodic() {return false;} */

/* 	BlockLabHITCloud_special(): BlockLab<BlockType,allocator>(){} */

/* 	void _apply_bc(const BlockInfo& info, const Real t=0) */
/* 	{ */
/*         BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock); */

/*         ElementTypeBlock b; */
/*         b.clear(); */
/*         b.rho = 1000; */
/*         b.u = b.v = b.w = 0; */
/*         b.G = 1./(6.12-1); */
/*         b.P = 3.43e3*b.G*6.12; */
/*         b.energy = 400.0; */

/*         if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(b); */
/*         if (info.index[0]==this->NX-1)  bc.template applyBC_dirichlet<0,1>(b); */
/*         if (info.index[1]==0)			bc.template applyBC_dirichlet<1,0>(b); */
/*         if (info.index[1]==this->NY-1)	bc.template applyBC_dirichlet<1,1>(b); */
/*         if (info.index[2]==0)			bc.template applyBC_dirichlet<2,0>(b); */
/*         if (info.index[2]==this->NZ-1)	bc.template applyBC_dirichlet<2,1>(b); */
/*     } */
/* }; */
/* typedef BlockLabHITCloud_special<Block_t, std::allocator> Lab; */

/* // cloud special dirichlet */
/* #include "Test_Cloud.h" */
/* template<typename BlockType, template<typename X> class allocator=std::allocator> */
/* class BlockLabCloudLaplace_special: public BlockLab<BlockType,allocator> */
/* { */
/* 	typedef typename BlockType::ElementType ElementTypeBlock; */

/* public: */

/*     virtual inline std::string name() const { return "BlockLabCloudLaplace_special"; } */
/* 	bool is_xperiodic() {return false;} */
/* 	bool is_yperiodic() {return false;} */
/* 	bool is_zperiodic() {return false;} */

/* 	BlockLabCloudLaplace_special(): BlockLab<BlockType,allocator>(){} */

/* 	void _apply_bc(const BlockInfo& info, const Real t=0) */
/* 	{ */
/*         BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock); */

/*         ElementTypeBlock b; */
/*         b.clear(); */
/*         b.rho = 1000; */
/*         b.u = b.v = b.w = 0; */
/*         b.G = 1./(6.12-1); */
/*         b.P = 3.430e3*b.G*6.12; */
/*         b.energy = 400.0; */

/*         if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(b); */
/*         if (info.index[0]==this->NX-1)  bc.template applyBC_dirichlet<0,1>(b); */
/*         if (info.index[1]==0)			bc.template applyBC_dirichlet<1,0>(b); */
/*         if (info.index[1]==this->NY-1)	bc.template applyBC_dirichlet<1,1>(b); */
/*         if (info.index[2]==0)			bc.template applyBC_dirichlet<2,0>(b); */
/*         if (info.index[2]==this->NZ-1)	bc.template applyBC_dirichlet<2,1>(b); */
/*     } */
/* }; */
/* typedef BlockLabCloudLaplace_special<Block_t, std::allocator> Lab; */


/* // hit cloud special dirichlet inflow */
/* #include "Test_Cloud.h" */
/* #include "Test_HIT.h" */
/* template<typename BlockType, template<typename X> class allocator=std::allocator> */
/* class BlockLabHITCloudInflow: public BlockLab<BlockType,allocator> */
/* { */
/* 	typedef typename BlockType::ElementType ElementTypeBlock; */

/* public: */

/*     virtual inline std::string name() const { return "BlockLabHITCloudInflow"; } */
/* 	bool is_xperiodic() {return true;} */
/* 	bool is_yperiodic() {return true;} */
/* 	bool is_zperiodic() {return false;} */

/* 	BlockLabHITCloudInflow(): BlockLab<BlockType,allocator>(){} */

/* 	void _apply_bc(const BlockInfo& info, const Real t=0) */
/* 	{ */
/*         BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock); */

/*         ElementTypeBlock b; */
/*         b.clear(); */
/*         b.G = 1./(6.12-1); */
/*         b.P = 3.430e3*b.G*6.12; */
/*         const double rho = HITData::rho1; */
/*         const double ru = rho * HITData::u0; */
/*         const double rv = rho * HITData::v0; */
/*         const double rw = rho * HITData::w0; */
/*         const double energy = b.G*HITData::p1 + b.P + 0.5 * (ru*ru + rv*rv + rw*rw)/rho; */
/*         b.rho = rho; */
/*         b.u = ru; */
/*         b.v = rv; */
/*         b.w = rw; */
/*         b.energy = energy; */

/*         if (info.index[2]==0)           bc.template applyBC_dirichlet<2,0>(b); */
/*         if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>(); */
/*     } */
/* }; */

/* template<typename BlockType, template<typename X> class allocator=std::allocator> */
/* class BlockLabHITCloudInflowPressureOut: public BlockLab<BlockType,allocator> */
/* { */
/* 	typedef typename BlockType::ElementType ElementTypeBlock; */

/* public: */

/*     virtual inline std::string name() const { return "BlockLabHITCloudInflowPressureOut"; } */
/* 	bool is_xperiodic() {return true;} */
/* 	bool is_yperiodic() {return true;} */
/* 	bool is_zperiodic() {return false;} */

/* 	BlockLabHITCloudInflowPressureOut(): BlockLab<BlockType,allocator>(){} */

/* 	void _apply_bc(const BlockInfo& info, const Real t=0) */
/* 	{ */
/*         BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock); */

/*         ElementTypeBlock b; */
/*         b.clear(); */
/*         b.G = 1./(6.12-1); */
/*         b.P = 3.430e3*b.G*6.12; */
/*         const double rho = HITData::rho1; */
/*         const double ru = rho * HITData::u0; */
/*         const double rv = rho * HITData::v0; */
/*         const double rw = rho * HITData::w0; */
/*         const double energy = b.G*HITData::p1 + b.P + 0.5 * (ru*ru + rv*rv + rw*rw)/rho; */
/*         b.rho = rho; */
/*         b.u = ru; */
/*         b.v = rv; */
/*         b.w = rw; */
/*         b.energy = energy; */

/*         if (info.index[2]==0)           bc.template applyBC_dirichlet_inflow<2,0>(b); */
/*         if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing_constPressure<2,1>(HITData::p1); */
/*     } */
/* }; */
/* typedef BlockLabHITCloudInflow<Block_t, std::allocator> Lab; */
/* /1* typedef BlockLabHITCloudInflowPressureOut<Block_t, std::allocator> Lab; *1/ */

// Kelvin-Helmholtz instability
/* #include "Test_KHinstability.h" */
/* typedef BlockLabKHinstability<Block_t, std::allocator> Lab; */
