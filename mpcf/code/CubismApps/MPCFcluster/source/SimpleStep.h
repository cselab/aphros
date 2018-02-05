#pragma once

#include <BlockLabMPI.h>
#include "BlockProcessor_MPI.h"
#include "Types.h"
#include "Tests.h"
#include <StencilInfo.h>
#include <Diffusion_CPP.h>

typedef BlockLabMPI<Lab> LabMPI;
typedef GridMPI< Grid_t > GridMPI_t;

namespace SimpleData
{
  int step_id = 0;

  // Diffusion is an adapter to call Kernel
  // Common usage:
  // - TKernel (process<TLab, TKernel, TGrid>() in BlockProcessor_MPI.h
  // - Processing (GridMPI<TGrid>::sync<Processing>() in GridMPI.h
  template < typename Kernel, typename Lab >
  struct Diffusion
  {
    StencilInfo stencil;
    Real dtinvh;
    int stencil_start[3];
    int stencil_end[3];

    StencilInfo getStencil() {
      return StencilInfo(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6);
    }

    Diffusion(const Real _dtinvh)
      : dtinvh(_dtinvh), stencil(getStencil())
    {
      std::cerr << "Diffusion::constructor(dtinvh)" << std::endl;
      stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
      stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }

    Diffusion(const Diffusion& c)
      : dtinvh(c.dtinvh), stencil(c.stencil)
    {
      std::cerr << "Diffusion::copy constructor" << std::endl;
      stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
      stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }

    inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
    {
      // At this point the only information available is:
      // - Lab, access to FluidElements via operator()(x,y,z)
      // - BlockInfo (defined in BlockInfo.h) containing index, origin, spacing
      // - Block_t=FluidBlock, 3D array with fields data and tmp
      // - dtinvh passed at construction
      //
      // Further parameters (e.g. viscosity) would be passed at construction.
      //
      // Affects: implementation of process<>()
      // Depends: interface of Kernel constructor and operator()

      // Create new instance of Kernel,
      // one instance per rank-step-block
      Kernel kernel(dtinvh);

      const Real * const srcfirst = &lab(-1,-1,-1).alpha2;
      const int labSizeRow = lab.template getActualSize<0>();
      const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
      Real * const destfirst =  &o.tmp[0][0][0][0];
      // Call kernel evaluation 
      kernel.compute(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
          destfirst, Block_t::gptfloats, 
          Block_t::sizeX, Block_t::sizeX*Block_t::sizeY);
    }
  };

}

template<typename TGrid>
class SimpleStep 
{
  TGrid & grid;

 public:
  ~SimpleStep() {}

  SimpleStep(TGrid& grid)
    : grid(grid)
  {}

  Real operator()(const Real dt, const Real current_time=0.0) {
    // At this point the only information available is:
    // - dt, t 
    // - TGrid=GridMPI=GridMPI<Grid<FluidBlock>>, 3D array for blocks
    //
    // Affects: implementation of Test_Simple
    // Depends: interface of Diffusion constructor and process<>()

    MPI_Comm comm = grid.getCartComm();
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    const bool isroot = (0 == myrank);

    std::vector<BlockInfo> vInfo = grid.getBlocksInfo();

    const Real h = vInfo[0].h_gridpoint; // uniform grid size
    const Real dtinvh = dt / h;

    // diffusion
    // Create a new instance of Diffusion, one intance per rank-step
    SimpleData::Diffusion<Diffusion_CPP, Lab> diffusion(dtinvh);
    process< LabMPI >(diffusion, (GridMPI_t&)grid, current_time, 0);

    // update
    for (size_t r = 0; r < vInfo.size(); ++r) {
      BlockInfo& bi = vInfo[r];
      Block_t& b = *(Block_t*)bi.ptrBlock;
      const Real* src = &b.tmp[0][0][0][0];
      Real* dst = &b.data[0][0][0].alpha2;
      int fe = Block_t::gptfloats;  // floats per element
      int n = std::pow(_BLOCKSIZE_, 3) * fe; // floats per block
      for (int i = 0; i < n; i += fe) {
        for (int k = 0; k < fe; ++k) {
          dst[i+k] += src[i+k];
          assert(isfinite(src[i+k]));
          assert(isfinite(dst[i+k]));
        }
      }
    }

    SimpleData::step_id++;

    return dt;
  }
};
