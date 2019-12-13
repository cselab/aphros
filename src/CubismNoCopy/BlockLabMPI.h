/*
 *  BlockLabMPI.h
 *  Cubism (allocation free)
 *
 *  Created by Fabian Wermelinger on 04/11/19.
 *  Copyright 2019 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#ifndef BLOCKLABMPI_H_6UGGRS5A
#define BLOCKLABMPI_H_6UGGRS5A

#include <stdio.h>

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-04] Get rid of this. Careful with down
// casts (due to inheritance)
#include "GridMPI.h"

template <typename MyBlockLab>
class BlockLabMPI : public MyBlockLab {
  using BlockType = typename MyBlockLab::BlockType;
  using Synch = SynchronizerMPI<BlockType>;
  const Synch* refSynchronizerMPI;

 protected:
  int mypeindex[3], pesize[3], mybpd[3];
  int gLastX, gLastY, gLastZ;

 public:
  template <typename TGrid>
  void prepare(GridMPI<TGrid>& grid, const Synch& SynchronizerMPI) {
    refSynchronizerMPI = &SynchronizerMPI;
    refSynchronizerMPI->getpedata(mypeindex, pesize, mybpd);
    StencilInfo stencil = refSynchronizerMPI->getstencil();
    assert(stencil.isvalid());
    MyBlockLab::prepare(
        grid, stencil.sx, stencil.ex, stencil.sy, stencil.ey, stencil.sz,
        stencil.ez, stencil.tensorial);
    gLastX = grid.getBlocksPerDimension(0) - 1;
    gLastY = grid.getBlocksPerDimension(1) - 1;
    gLastZ = grid.getBlocksPerDimension(2) - 1;
  }

  template <typename BlockType>
  void load(BlockType& target, std::vector<BlockType>& vfields) {
    // global block index
    const int index[3] = {target.index[0], target.index[1], target.index[2]};

    MyBlockLab::load(target, vfields);

    assert(refSynchronizerMPI != NULL);

    const int xorigin = mypeindex[0] * mybpd[0];
    const int yorigin = mypeindex[1] * mybpd[1];
    const int zorigin = mypeindex[2] * mybpd[2];

    const bool xskin =
        (index[0] == xorigin || index[0] == xorigin + mybpd[0] - 1);
    const bool yskin =
        (index[1] == yorigin || index[1] == yorigin + mybpd[1] - 1);
    const bool zskin =
        (index[2] == zorigin || index[2] == zorigin + mybpd[2] - 1);

    // const bool xboundary = index[0]==0 || index[0]==gLastX;
    // const bool yboundary = index[1]==0 || index[1]==gLastY;
    // const bool zboundary = index[2]==0 || index[2]==gLastZ;

    // const bool any_periodic = this->is_xperiodic() || this->is_yperiodic() ||
    // this->is_zperiodic(); const bool any_boundary = xboundary || yboundary ||
    // zboundary;
    if ((xskin || yskin ||
         zskin)) // && (!any_boundary || any_boundary && any_periodic))
    {
      const bool xperiodic = this->is_xperiodic();
      const bool yperiodic = this->is_yperiodic();
      const bool zperiodic = this->is_zperiodic();

      const int rsx =
          (!xperiodic && index[0] == 0) ? 0 : this->m_stencilStart[0];
      const int rex = (!xperiodic && index[0] == gLastX)
                          ? BlockType::sizeX
                          : (BlockType::sizeX + this->m_stencilEnd[0] - 1);
      const int rsy =
          (!yperiodic && index[1] == 0) ? 0 : this->m_stencilStart[1];
      const int rey = (!yperiodic && index[1] == gLastY)
                          ? BlockType::sizeY
                          : (BlockType::sizeY + this->m_stencilEnd[1] - 1);
      const int rsz =
          (!zperiodic && index[2] == 0) ? 0 : this->m_stencilStart[2];
      const int rez = (!zperiodic && index[2] == gLastZ)
                          ? BlockType::sizeZ
                          : (BlockType::sizeZ + this->m_stencilEnd[2] - 1);

      refSynchronizerMPI->fetch_soa(target, rsx, rex, rsy, rey, rsz, rez);
    }

    // duplicate slices for 2D case
    if (1 == BlockType::sizeZ) {
      this->fix2D(target);
    }
  }
};

#endif /* BLOCKLABMPI_H_6UGGRS5A */
