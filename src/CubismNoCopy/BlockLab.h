/*
 *  BlockLabMPI.h
 *  Cubism (allocation free)
 *
 *  Created by Fabian Wermelinger on 04/11/19.
 *  Copyright 2019 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#ifndef BLOCKLAB_H_PRFDFEI9
#define BLOCKLAB_H_PRFDFEI9

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include <cstring>
#include <string>
#include <vector>

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-04] Get rid of this include!
#include "Grid.h"

#ifdef __bgq__
#include <builtins.h>
#define memcpy2(a, b, c) __bcopy((b), (a), (c))
#else
#define memcpy2(a, b, c) std::memcpy((a), (b), (c))
#endif

/**
 * Working copy of Block + Ghosts.
 * Data of original block is copied (!) here. So when changing something in
 * the lab we are not changing the original data.
 * Requirements:
 * - Concepts::BlockLabElementType<ElementTypeT>
 * - Concepts::Castable<typename BlockType::ElementType, ElementTypeT>
 */
template <typename TBlock>
class BlockLab {
  // check concepts
  //   CONCEPT_CHECK(Concepts::BlockLabElementType<ElementTypeT>);
  //   CONCEPT_CHECK(Concepts::Castable<typename BlockType::ElementType,
  //   ElementTypeT>);

 public:
  typedef typename TBlock::ElementType ElementType;

 protected:
  typedef TBlock BlockType;
  typedef typename BlockType::ElementType ElementTypeBlock;
  enum eBlockLab_State {
    eMRAGBlockLab_Prepared,
    eMRAGBlockLab_Loaded,
    eMRAGBlockLab_Uninitialized
  };

  eBlockLab_State m_state;

  int m_stencilStart[3], m_stencilEnd[3];
  int NX, NY, NZ;

  bool istensorial;

  const Grid<BlockType>* m_refGrid;

  void fix2D(TBlock& target) {
    using TReal = typename TBlock::ElementType;
    const int zslice = target.n_comp * TBlock::stridex * TBlock::stridey;
    const size_t bytes = zslice * sizeof(TReal);
    const char* src = (const char*)(target.data_halo + TBlock::halo * zslice);
    TReal* const dst_lo = target.data_halo;
    TReal* const dst_hi = target.data_halo + (TBlock::halo + 1) * zslice;
    for (int i = 0; i < TBlock::halo; ++i) {
      memcpy2((char*)(dst_lo + i * zslice), src, bytes);
      memcpy2((char*)(dst_hi + i * zslice), src, bytes);
    }
  }

 public:
  BlockLab() : m_state(eMRAGBlockLab_Uninitialized), m_refGrid(NULL) {
    m_stencilStart[0] = m_stencilStart[1] = m_stencilStart[2] = 0;
    m_stencilEnd[0] = m_stencilEnd[1] = m_stencilEnd[2] = 0;
  }

  BlockLab(const BlockLab& c) = delete;
  BlockLab& operator=(const BlockLab& c) = delete;

  virtual std::string name() const {
    return "BlockLab";
  }
  virtual bool is_xperiodic() {
    return true;
  }
  virtual bool is_yperiodic() {
    return true;
  }
  virtual bool is_zperiodic() {
    return true;
  }

  virtual ~BlockLab() {}

  void prepare(
      Grid<BlockType>& grid, int startX, int endX, int startY, int endY,
      int startZ, int endZ, const bool istensorial) {
    const int ss[3] = {startX, startY, startZ};
    const int se[3] = {endX, endY, endZ};
    prepare(grid, ss, se, istensorial);
  }

  /**
   * Prepare the extended block.
   * @param collection    Collection of blocks in the grid (e.g. result of
   * Grid::getBlockCollection()).
   * @param boundaryInfo  Info on the boundaries of the grid (e.g. result of
   * Grid::getBoundaryInfo()).
   * @param stencil_start Maximal stencil used for computations at lower
   * boundary. Defines how many ghosts we will get in extended block.
   * @param stencil_end   Maximal stencil used for computations at lower
   * boundary. Defines how many ghosts we will get in extended block.
   */

  void prepare(
      Grid<BlockType>& grid, const int stencil_start[3],
      const int stencil_end[3], const bool istensorial) {
    NX = grid.getBlocksPerDimension(0);
    NY = grid.getBlocksPerDimension(1);
    NZ = grid.getBlocksPerDimension(2);

    this->istensorial = istensorial;

    m_refGrid = &grid;

    assert(stencil_start[0] >= -BlockType::sizeX);
    assert(stencil_start[1] >= -BlockType::sizeY);
    assert(stencil_start[2] >= -BlockType::sizeZ);
    assert(stencil_end[0] < BlockType::sizeX * 2);
    assert(stencil_end[1] < BlockType::sizeY * 2);
    // assert(stencil_end[2] < BlockType::sizeZ*2);

    m_stencilStart[0] = stencil_start[0];
    m_stencilStart[1] = stencil_start[1];
    m_stencilStart[2] = stencil_start[2];

    m_stencilEnd[0] = stencil_end[0];
    m_stencilEnd[1] = stencil_end[1];
    m_stencilEnd[2] = stencil_end[2];

    assert(m_stencilStart[0] <= m_stencilEnd[0]);
    assert(m_stencilStart[1] <= m_stencilEnd[1]);
    assert(m_stencilStart[2] <= m_stencilEnd[2]);

    m_state = eMRAGBlockLab_Prepared;
  }

  /**
   * Load a block (incl. ghosts for it).
   * This is not called internally but by the BlockProcessing-class. Hence a new
   * version of BlockLab, can just overwrite it and through template-passing to
   * BlockProcessing, the right version will be called.
   * @param info  Reference to info of block to be loaded.
   */

  void load(TBlock& target, std::vector<TBlock>& vfields) {
    // global block index
    const int index[3] = {target.index[0], target.index[1], target.index[2]};

    const Grid<BlockType>& grid = *m_refGrid; // holds neighbor information only

    // 0. couple of checks
    // 1. put the ghosts into the cache
    // 2. check if 2D block

    // 0.
    assert(
        m_state == eMRAGBlockLab_Prepared || m_state == eMRAGBlockLab_Loaded);

    const int nX = BlockType::sizeX;
    const int nY = BlockType::sizeY;
    const int nZ = BlockType::sizeZ;

    // 1.
    {
      const bool xperiodic = is_xperiodic();
      const bool yperiodic = is_yperiodic();
      const bool zperiodic = is_zperiodic();

      const bool xskin =
          index[0] == 0 || index[0] == (int)grid.getBlocksPerDimension(0) - 1;
      const bool yskin =
          index[1] == 0 || index[1] == (int)grid.getBlocksPerDimension(1) - 1;
      const bool zskin =
          index[2] == 0 || index[2] == (int)grid.getBlocksPerDimension(2) - 1;

      const int xskip = index[0] == 0 ? -1 : 1;
      const int yskip = index[1] == 0 ? -1 : 1;
      const int zskip = index[2] == 0 ? -1 : 1;

      for (int icode = 0; icode < 27; icode++) {
        if (icode == 1 * 1 + 3 * 1 + 9 * 1) continue;

        const int code[3] = {
            icode % 3 - 1, (icode / 3) % 3 - 1, (icode / 9) % 3 - 1};

        if (!xperiodic && code[0] == xskip && xskin) continue;
        if (!yperiodic && code[1] == yskip && yskin) continue;
        if (!zperiodic && code[2] == zskip && zskin) continue;

        if (!istensorial && abs(code[0]) + abs(code[1]) + abs(code[2]) > 1)
          continue;

        const int s[3] = {
            code[0] < 1 ? (code[0] < 0 ? m_stencilStart[0] : 0) : nX,
            code[1] < 1 ? (code[1] < 0 ? m_stencilStart[1] : 0) : nY,
            code[2] < 1 ? (code[2] < 0 ? m_stencilStart[2] : 0) : nZ};

        const int e[3] = {
            code[0] < 1 ? (code[0] < 0 ? 0 : nX) : nX + m_stencilEnd[0] - 1,
            code[1] < 1 ? (code[1] < 0 ? 0 : nY) : nY + m_stencilEnd[1] - 1,
            code[2] < 1 ? (code[2] < 0 ? 0 : nZ) : nZ + m_stencilEnd[2] - 1};

        BlockType& b = vfields[grid(
            index[0] + code[0], index[1] + code[1], index[2] + code[2])];
#if 1
        const int m_vSize0 = TBlock::stridex;
        const int m_nElemsPerSlice = TBlock::stridex * TBlock::stridey;

        const int my_ix = s[0] - m_stencilStart[0];

        // printf("iy : %d %d\n", s[1], e[1]);
        const int bytes = (e[0] - s[0]) * target.n_comp * sizeof(ElementType);
        for (int iz = s[2]; iz < e[2]; iz++) {
          const int my_izx =
              (iz - m_stencilStart[2]) * m_nElemsPerSlice + my_ix;
#if 0
					for(int iy=s[1]; iy<e[1]; iy++)
					{
#if 1 // ...
      // char * ptrDest = (char*)&m_cacheBlock->Access(s[0]-m_stencilStart[0],
      // iy-m_stencilStart[1], iz-m_stencilStart[2]);
						char * ptrDest = (char*)&m_cacheBlock->LinAccess(my_izx + (iy-m_stencilStart[1])*m_vSize0);

						const char * ptrSrc = (const char*)&b(s[0] - code[0]*BlockType::sizeX, iy - code[1]*BlockType::sizeY, iz - code[2]*BlockType::sizeZ);
						memcpy2((char *)ptrDest, (char *)ptrSrc, bytes);
#else
						for(int ix=s[0]; ix<e[0]; ix++)
							target.LinAccess(ix-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]) =
							(ElementType)b(ix - code[0]*BlockType::sizeX, iy - code[1]*BlockType::sizeY, iz - code[2]*BlockType::sizeZ);
#endif
					}
#else
          if ((e[1] - s[1]) % 4 != 0) {
            for (int iy = s[1]; iy < e[1]; iy++) {
              char* ptrDest = (char*)&target.LinAccess(
                  my_izx + (iy - m_stencilStart[1]) * m_vSize0);

              // assert(b(s[0] - code[0]*BlockType::sizeX, iy -
              // code[1]*BlockType::sizeY, iz - code[2]*BlockType::sizeZ) > 0);

              // XXX: [fabianw@mavt.ethz.ch; 2019-11-11] This is
              // dangerous if b() returns something complex that
              // is not POD!
              const char* ptrSrc = (const char*)&b(
                  s[0] - code[0] * BlockType::sizeX,
                  iy - code[1] * BlockType::sizeY,
                  iz - code[2] * BlockType::sizeZ);
              memcpy2((char*)ptrDest, (char*)ptrSrc, bytes);
            }
          } else {
            for (int iy = s[1]; iy < e[1]; iy += 4) {
              char* ptrDest0 = (char*)&target.LinAccess(
                  my_izx + (iy + 0 - m_stencilStart[1]) * m_vSize0);
              char* ptrDest1 = (char*)&target.LinAccess(
                  my_izx + (iy + 1 - m_stencilStart[1]) * m_vSize0);
              char* ptrDest2 = (char*)&target.LinAccess(
                  my_izx + (iy + 2 - m_stencilStart[1]) * m_vSize0);
              char* ptrDest3 = (char*)&target.LinAccess(
                  my_izx + (iy + 3 - m_stencilStart[1]) * m_vSize0);

              // XXX: [fabianw@mavt.ethz.ch; 2019-11-11] This is
              // dangerous if b() returns something complex that
              // is not POD!
              const char* ptrSrc0 = (const char*)&b(
                  s[0] - code[0] * BlockType::sizeX,
                  iy + 0 - code[1] * BlockType::sizeY,
                  iz - code[2] * BlockType::sizeZ);
              const char* ptrSrc1 = (const char*)&b(
                  s[0] - code[0] * BlockType::sizeX,
                  iy + 1 - code[1] * BlockType::sizeY,
                  iz - code[2] * BlockType::sizeZ);
              const char* ptrSrc2 = (const char*)&b(
                  s[0] - code[0] * BlockType::sizeX,
                  iy + 2 - code[1] * BlockType::sizeY,
                  iz - code[2] * BlockType::sizeZ);
              const char* ptrSrc3 = (const char*)&b(
                  s[0] - code[0] * BlockType::sizeX,
                  iy + 3 - code[1] * BlockType::sizeY,
                  iz - code[2] * BlockType::sizeZ);

              memcpy2((char*)ptrDest0, (char*)ptrSrc0, bytes);
              memcpy2((char*)ptrDest1, (char*)ptrSrc1, bytes);
              memcpy2((char*)ptrDest2, (char*)ptrSrc2, bytes);
              memcpy2((char*)ptrDest3, (char*)ptrSrc3, bytes);
            }
          }
#endif
        }
#else
        const int off_x = -code[0] * nX + m_stencilStart[0];
        const int off_y = -code[1] * nY + m_stencilStart[1];
        const int off_z = -code[2] * nZ + m_stencilStart[2];

        const int nbytes = (e[0] - s[0]) * target.n_comp * sizeof(ElementType);
#if 1
        const int _iz0 = s[2] - m_stencilStart[2];
        const int _iz1 = e[2] - m_stencilStart[2];
        const int _iy0 = s[1] - m_stencilStart[1];
        const int _iy1 = e[1] - m_stencilStart[1];

        for (int iz = _iz0; iz < _iz1; iz++)
          for (int iy = _iy0; iy < _iy1; iy++)

#else
        for (int iz = s[2] - m_stencilStart[2]; iz < e[2] - m_stencilStart[2];
             iz++)
          for (int iy = s[1] - m_stencilStart[1]; iy < e[1] - m_stencilStart[1];
               iy++)
#endif
          {
#if 0
						char * ptrDest = (char*)&m_cacheBlock->Access(s[0]-m_stencilStart[0], iy, iz);
						const char * ptrSrc = (const char*)&b(0 + off_x, iy + off_y, iz + off_z);
						memcpy2(ptrDest, ptrSrc, nbytes);
#else
            for (int ix = s[0] - m_stencilStart[0];
                 ix < e[0] - m_stencilStart[0]; ix++)
              m_cacheBlock->Access(ix, iy, iz) =
                  (ElementType)b(ix + off_x, iy + off_y, iz + off_z);
#endif
          }
#endif
      }

      m_state = eMRAGBlockLab_Loaded;
    }

    // 2.
    if (1 == TBlock::sizeZ) {
      fix2D(target);
    }
  }
};

#endif /* BLOCKLAB_H_PRFDFEI9 */
