/*
 *  Grid.h
 *  Cubism (allocation free)
 *
 *  Created by Fabian Wermelinger on 04/11/19.
 *  Copyright 2019 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#ifndef GRID_H_I6LVUQY7
#define GRID_H_I6LVUQY7

#include <cassert>
#include <vector>

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-08] get rid of this include (template)
#include "BlockInfo.h"

template <typename B_>
class Grid {
 public:
  using Block = B_;
  using BlockType = Block;

 private:
  // FIXME: [fabianw@mavt.ethz.ch; 2019-11-04] BlockInfo type could be
  // replaced with templated type.
  std::vector<BlockInfo> m_vInfo;

 protected:
  const size_t NX, NY, NZ, N;
  const double maxextent;

  size_t _encode(const size_t ix, const size_t iy, const size_t iz) const {
    assert(ix < NX);
    assert(iy < NY);
    assert(iz < NZ);

    return ix + NX * (iy + NY * iz);
  }

 public:
  Grid(
      const size_t NX_, const size_t NY_ = 1, const size_t NZ_ = 1,
      const double maxextent_ = 1)
      : NX(NX_), NY(NY_), NZ(NZ_), N(NX_ * NY_ * NZ_), maxextent(maxextent_) {
    // cells per block
    const size_t bs[3] = {Block::sizeX, Block::sizeY, Block::sizeZ};
    // cells total
    const size_t nn[3] = {NX * bs[0], NY * bs[1], NZ * bs[2]};
    // cell size (h_gridpoint from BlockInfo)
    const double hc = (maxextent / std::max(nn[0], std::max(nn[1], nn[2])));
    // block extent
    const double h = std::max(bs[0], std::max(bs[1], bs[2])) * hc;

    for (size_t iz = 0; iz < NZ; iz++)
      for (size_t iy = 0; iy < NY; iy++)
        for (size_t ix = 0; ix < NX; ix++) {
          const long long blockID = _encode(ix, iy, iz);
          const int idx[3] = {(int)ix, (int)iy, (int)iz};
          const double origin[3] = {ix * h, iy * h, iz * h};
          m_vInfo.push_back(BlockInfo(blockID, idx, origin, h, hc));
        }
  }

  virtual ~Grid() {}

  virtual size_t getBlocksPerDimension(int idim) const {
    assert(idim >= 0 && idim < 3);

    switch (idim) {
      case 0:
        return NX;
      case 1:
        return NY;
      case 2:
        return NZ;
      default:
        abort();
        return 0;
    }
  }

  virtual size_t operator()(int ix, int iy = 0, int iz = 0) const {
    return _encode((ix + NX) % NX, (iy + NY) % NY, (iz + NZ) % NZ);
  }

  virtual std::vector<BlockInfo>& getBlocksInfo() {
    return m_vInfo;
  }

  virtual const std::vector<BlockInfo>& getBlocksInfo() const {
    return m_vInfo;
  }
};

#endif /* GRID_H_I6LVUQY7 */
