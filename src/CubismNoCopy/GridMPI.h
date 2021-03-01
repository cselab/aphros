/*
 *  GridMPI.h
 *  Cubism (allocation free)
 *
 *  Created by Fabian Wermelinger on 04/11/19.
 *  Copyright 2019 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#ifndef GRIDMPI_H_J17BPEDT
#define GRIDMPI_H_J17BPEDT

#include <mpi.h>
#include <map>
#include <string>
#include <vector>

#include "StencilInfo.h"
#include "SynchronizerMPI.h"

template <typename TGrid>
class GridMPI : public TGrid {
 public:
  using Block = typename TGrid::Block;
  using BlockType = typename TGrid::BlockType;
  using Synch = SynchronizerMPI<Block>;
  friend Synch;

 private:
  size_t timestamp;

 protected:
  int myrank, mypeindex[3], pesize[3];
  int periodic[3];
  int mybpd[3], myblockstotalsize, blocksize[3];

  std::vector<BlockInfo> cached_blockinfo;

  // acts like cache
  std::map<StencilInfo, Synch*> SynchronizerMPIs;

  MPI_Comm worldcomm;
  MPI_Comm cartcomm;

 public:
  GridMPI(
      const int npeX, const int npeY, const int npeZ, const int nX,
      const int nY = 1, const int nZ = 1, const double maxextent = 1,
      const MPI_Comm comm = MPI_COMM_WORLD)
      : TGrid(nX, nY, nZ, maxextent), timestamp(0), worldcomm(comm) {
    blocksize[0] = Block::sizeX;
    blocksize[1] = Block::sizeY;
    blocksize[2] = Block::sizeZ;

    mybpd[0] = nX;
    mybpd[1] = nY;
    mybpd[2] = nZ;
    myblockstotalsize = nX * nY * nZ;

    periodic[0] = true;
    periodic[1] = true;
    periodic[2] = true;

    pesize[0] = npeX;
    pesize[1] = npeY;
    pesize[2] = npeZ;

    int world_size;
    MPI_Comm_size(worldcomm, &world_size);
    assert(npeX * npeY * npeZ == world_size);

    MPI_Cart_create(worldcomm, 3, pesize, periodic, true, &cartcomm);
    MPI_Comm_rank(cartcomm, &myrank);
    MPI_Cart_coords(cartcomm, myrank, 3, mypeindex);

    const std::vector<BlockInfo> vInfo = TGrid::getBlocksInfo();

    // cells total
    const int nn[3] = {
        mybpd[0] * blocksize[0] * pesize[0],
        mybpd[1] * blocksize[1] * pesize[1],
        mybpd[2] * blocksize[2] * pesize[2]};
    // cell size (h_gridpoint from BlockInfo)
    const double hc = (maxextent / std::max(nn[0], std::max(nn[1], nn[2])));
    // block extent
    const double h =
        std::max(blocksize[0], std::max(blocksize[1], blocksize[2])) * hc;

    for (size_t i = 0; i < vInfo.size(); ++i) {
      BlockInfo info = vInfo[i];
      info.h = h;
      info.h_gridpoint = hc;
      for (int j = 0; j < 3; ++j) {
        info.index[j] += mypeindex[j] * mybpd[j];
        info.origin[j] = info.index[j] * h;
      }
      cached_blockinfo.push_back(info);
    }
  }

  ~GridMPI() {
    for (typename std::map<StencilInfo, Synch*>::const_iterator it =
             SynchronizerMPIs.begin();
         it != SynchronizerMPIs.end(); ++it)
      delete it->second;

    SynchronizerMPIs.clear();
    MPI_Comm_free(&cartcomm);
  }

  std::vector<BlockInfo>& getBlocksInfo() override {
    return cached_blockinfo;
  }

  const std::vector<BlockInfo>& getBlocksInfo() const override {
    return cached_blockinfo;
  }

  std::vector<BlockInfo>& getResidentBlocksInfo() {
    return TGrid::getBlocksInfo();
  }

  const std::vector<BlockInfo>& getResidentBlocksInfo() const {
    return TGrid::getBlocksInfo();
  }

  // for a given kernel p (e.g. Processing=Diffusion)
  // based on StencilInfo (bounding box and selected components)
  // returns a SynchronizerMPI (new or existing)
  // and performs communication
  // (e.g. two kernels with identical stencils and
  // same elected components would share a common SynchronizerMPI)
  template <typename TView>
  Synch& sync(
      std::vector<std::vector<TView>>& fields, const int nhalo_start[3],
      const int nhalo_end[3], const bool is_tensorial = true,
      const bool compress = true, const bool hist = true) {
    int nhalo_sz = nhalo_start[2];
    int nhalo_ez = nhalo_end[2];
    if (1 == TView::sizeZ) {
      // communication needed for tensorial
      nhalo_sz = 1;
      nhalo_ez = 1;
    }
    StencilInfo stencil(
        -nhalo_start[0], -nhalo_start[1], -nhalo_sz, nhalo_end[0] + 1,
        nhalo_end[1] + 1, nhalo_ez + 1, is_tensorial, 1, 0);
    int total_components = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
      const auto& f = fields[i];
      total_components += f[0].n_comp;
#if !defined(NDEBUG)
      for (size_t j = 1; j < f.size(); ++j) {
        assert(f[j - 1].n_comp == f[j].n_comp);
      }
#endif
    }
    stencil.selcomponents.clear();
    for (int i = 0; i < total_components; ++i) {
      stencil.selcomponents.push_back(i);
    }
    assert(stencil.isvalid());

    Synch* queryresult = NULL;

    typename std::map<StencilInfo, Synch*>::iterator itSynchronizerMPI =
        SynchronizerMPIs.find(stencil);

    if (itSynchronizerMPI == SynchronizerMPIs.end()) {
#if !defined(NDEBUG)
      printf(
          "Alloc sync: rank=%d; nfields=%lu; s=(%d,%d,%d); "
          "e=(%d,%d,%d); tensorial=%d\n",
          myrank, fields.size(), nhalo_start[0], nhalo_start[1], nhalo_start[2],
          nhalo_end[0], nhalo_end[1], nhalo_end[2], is_tensorial);
#endif
      queryresult = new Synch(
          fields, getBlocksInfo(), stencil, cartcomm, mybpd, blocksize,
          compress, hist);

      SynchronizerMPIs[stencil] = queryresult;
    } else {
      queryresult = itSynchronizerMPI->second;
    }

    // perform communication
    queryresult->sync(fields, timestamp);

    timestamp = (timestamp + 1) % 32768;

    return *queryresult;
  }

  Synch& getSynchronizerMPI(const StencilInfo stencil) {
    assert((SynchronizerMPIs.find(stencil) != SynchronizerMPIs.end()));

    return *SynchronizerMPIs.find(stencil)->second;
  }

  const std::map<StencilInfo, Synch*>& getSynchronizerMPI() const {
    return SynchronizerMPIs;
  }

  size_t getResidentBlocksPerDimension(int idim) const {
    assert(idim >= 0 && idim < 3);
    return (size_t)mybpd[idim];
  }

  size_t getBlocksPerDimension(int idim) const override {
    assert(idim >= 0 && idim < 3);
    return (size_t)mybpd[idim] * pesize[idim];
  }

  void peindex(int mypeindex[3]) const {
    for (int i = 0; i < 3; ++i)
      mypeindex[i] = this->mypeindex[i];
  }

  size_t getTimeStamp() const {
    return timestamp;
  }

  MPI_Comm getCartComm() const {
    return cartcomm;
  }

  MPI_Comm getWorldComm() const {
    return worldcomm;
  }
};

#endif /* GRIDMPI_H_J17BPEDT */
