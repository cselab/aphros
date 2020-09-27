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
  static constexpr const char* kDefaultName = "data";
  // Writes HDF file with dimensions (nz, ny, nx, 1).
  // fc: input field
  // path: output path
  // dname: dataset name
  static void Write(
      const FieldCell<typename M::Scal>& fc, std::string path, M& m,
      const std::string dname = kDefaultName);
  // Reads HDF file with one scalar field. Dimensions of the dataset
  // must match the global mesh size: (nz, ny, nx, 1).
  //
  // fc: output field
  // path: input path
  // dname: dataset name
  static void Read(
      FieldCell<typename M::Scal>& fc, std::string path, M& m,
      const std::string dname = kDefaultName);
  // Returns shape of a dataset in HDF file.
  // path: input path
  // dname: dataset name
  static std::vector<size_t> GetShape(
      std::string path, const std::string dname = kDefaultName);
  // Creates XMF file with uniform mesh and scalar data attribute linked to HDF.
  // xmfpath: output path
  // name: name of attribute
  // origin: mesh origin
  // spacing: mesh spacing (cell size)
  // dims: mesh size in cells (x,y,z)
  // hdfpath: path to existing hdf
  // dname: hdf dataset name
  static void WriteXmf(
      std::string xmfpath, std::string name,
      const std::array<double, 3>& origin, const std::array<double, 3>& spacing,
      const std::array<size_t, 3>& dims, std::string hdfpath,
      std::string dname);
  static void WriteXmf(
      std::string xmfpath, std::string name, std::string hdfpath, const M& m,
      std::string dname = kDefaultName);
};
