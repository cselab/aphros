// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <string>
#include "geom/mesh.h"

template <class M>
class Hdf {
 public:
  using Scal = typename M::Scal;
  using MIdx = typename M::MIdx;
  // Creates XMF file with uniform mesh and scalar data attribute
  // linked to HDF field with name 'data'.
  // xmfpath: path to output
  // name: name of attribute
  // origin,spacing,dims: mesh parameters
  // hdfpath: path to existing hdf
  // hdfdims: dimensions of hdf, either (nz,ny,nx,1) or (ny,nx,1)
  static void Write(
      const FieldCell<typename M::Scal>& fc, std::string path, M& m);
  // Creates XMF file with uniform mesh and scalar data attribute
  // linked to HDF field with name 'data'.
  // xmfpath: path to output
  // name: name of attribute
  // origin,spacing,dims: mesh parameters
  // hdfpath: path to existing hdf
  // hdfdims: dimensions of hdf, either (nz,ny,nx,1) or (ny,nx,1)
  static void WriteXmf(
      std::string xmfpath, std::string name,
      const std::array<double, 3>& origin, const std::array<double, 3>& spacing,
      const std::array<size_t, 3>& dims, std::string hdfpath);
  static void WriteXmf(
      std::string xmfpath, std::string name, std::string hdfpath, const M& m);
};
