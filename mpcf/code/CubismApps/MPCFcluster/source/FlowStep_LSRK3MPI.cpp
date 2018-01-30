/*
 *  FlowStep_LSRK3MPI.cpp
 *  MPCFcluster
 *
 *  Created by Babak Hejazi on 11/20/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <mpi.h>
#include "FlowStep_LSRK3MPI.h"

namespace LSRK3data {
    int verbosity;
    int step_id = 0;
}

template<typename TGrid>
template<typename Kdiff>
pair<Real, Real> FlowStep_LSRK3MPI<TGrid>::step(TGrid& grid, std::vector<BlockInfo>& vInfo, const Real a, const Real b, const Real dt, const Real current_time)
{
  const Real h = vInfo[0].h_gridpoint; // uniform grid size
  const Real dtinvh = dt / h;

  // diffusion
  LSRK3data::Diffusion<Kdiff, Lab> diffusion(dtinvh);
  process< LabMPI >(diffusion, (TGrid&)grid, current_time, 0);

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
  return pair<double, double>(1.,1.);
}

template<typename TGrid>
template<typename Kdiff>
inline void FlowStep_LSRK3MPI<TGrid>::_process_LSRK3(TGrid& grid, const Real dt, const Real current_time)
{
  std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
  step<Kdiff>(grid, vInfo, 1., 0., dt, current_time);
}

//Specialization
template<>
Real FlowStep_LSRK3MPI<GridMPI_t>::operator()(const Real dt, const Real current_time)
{
    MPI_Comm comm = grid.getCartComm();
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    const bool isroot = (0 == myrank);

    _process_LSRK3<Diffusion_CPP>(grid, dt, current_time);

    LSRK3data::step_id++;

    return dt;
}
