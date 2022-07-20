// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <string>
#include "geom/mesh.h"
#include "xmf.h"

namespace dump {

template <class M>
class Hdf {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  static constexpr const char* kDefaultName = "data";
  // Writes HDF file with dimensions (nz, ny, nx, 1).
  // fc: input field
  // path: output path
  // dname: dataset name
  template <class Field>
  static void Write(
      const Field& fc, std::string path, M& m,
      std::string dname = kDefaultName);
  // Writes blocks forming a multidimensional array to file.
  // path: path to HDF file
  // starts: index of start for each block
  // sizes: size of each block
  // data: data buffers for each block
  // global_size: total size of the array
  // type: type of entries to write
  // dname: dataset name
  // append: append existing file (true), create or truncate file (false)
  static void WriteBlocks(
      const std::string& path, const std::vector<MIdx>& starts,
      const std::vector<MIdx>& sizes,
      const std::vector<std::vector<Scal>>& data, MIdx global_size, Type type,
      std::string dname, const MpiWrapper& mpi, bool append);
  // Reads HDF file with one scalar field. Dimensions of the dataset
  // must match the global mesh size: (nz, ny, nx, 1).
  //
  // fc: output field
  // path: input path
  // dname: dataset name
  template <class Field>
  static void Read(
      Field& fc, std::string path, M& m, std::string dname = kDefaultName);
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
      std::string xmfpath, std::string name, Vect origin, Vect spacing,
      MIdx dims, std::string hdfpath, std::string dname);
  static void WriteXmf(
      std::string xmfpath, std::string name, std::string hdfpath, const M& m,
      std::string dname = kDefaultName);
};

} // namespace dump
