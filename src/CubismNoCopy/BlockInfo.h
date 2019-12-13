/*
 *  BlockInfo.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cstdlib>

struct BlockInfo {
  long long blockID;
  void* ptrBlock;
  bool special;
  int index[3];

  double origin[3];
  double h, h_gridpoint;
  double uniform_grid_spacing[3];
  double block_extent[3];

  double* ptr_grid_spacing[3];

  bool bUniform[3];

  template <typename T>
  inline void pos(T p[2], int ix, int iy) const {
    const int I[2] = {ix, iy};
    for (int j = 0; j < 2; ++j) {
      T delta = 0.0;
      if (bUniform[j])
        delta = uniform_grid_spacing[j] * (I[j] + 0.5);
      else {
        const double* const h = ptr_grid_spacing[j];
        for (int i = 0; i < I[j]; ++i)
          delta += h[i];
        delta += 0.5 * h[I[j]];
      }
      p[j] = origin[j] + delta;
    }
  }

  template <typename T>
  inline void pos(T p[3], int ix, int iy, int iz) const {
    const int I[3] = {ix, iy, iz};
    for (int j = 0; j < 3; ++j) {
      T delta = 0.0;
      if (bUniform[j])
        delta = uniform_grid_spacing[j] * (I[j] + 0.5);
      else {
        const double* const h = ptr_grid_spacing[j];
        for (int i = 0; i < I[j]; ++i)
          delta += h[i];
        delta += 0.5 * h[I[j]];
      }
      p[j] = origin[j] + delta;
    }
  }

  BlockInfo(
      long long ID, const int idx[3], const double pos[3], const double spacing,
      double h_gridpoint_, void* ptr = NULL, const bool _special = false)
      : blockID(ID), ptrBlock(ptr), special(_special) {
    h = spacing;
    h_gridpoint = h_gridpoint_;

    for (int i = 0; i < 3; ++i) {
      ptr_grid_spacing[i] = NULL;

      index[i] = idx[i];

      origin[i] = pos[i];

      uniform_grid_spacing[i] = h_gridpoint_;

      block_extent[i] = h;

      bUniform[i] = true;
    }
  }

  BlockInfo() : blockID(-1), ptrBlock(NULL) {}
};
