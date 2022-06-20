// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <istream>
#include <ostream>
#include <string>

#include "geom/mesh.h"
#include "xmf.h"

namespace dump {

template <class M>
class Raw {
 public:
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Xmf = dump::Xmf<Vect>;
  using Meta = typename Xmf::Meta;

  // Writes blocks forming a multidimensional array to file.
  // path: path to file
  // starts: index of start for each block
  // sizes: size of each block
  // data: data buffers for each block
  // global_size: total size of the array
  // type: type of entries to write
  // nompi: use a version without MPI even if supported
  template <class T>
  static void Write(
      const std::string& path, const std::vector<MIdx>& starts,
      const std::vector<MIdx>& sizes, const std::vector<std::vector<T>>& data,
      MIdx global_size, Type type, const MpiWrapper& mpi, bool nompi = false);
  template <class T>
  static void Write(
      const FieldCell<T>& fc, const Meta& meta, const std::string& path, M& m);
  // Writes field to raw-file and metadata to xmf-file.
  // fc: field to write
  // field_name: dataset name to put in metadata
  // raw_path: output path to file with extension `.raw`
  template <class T>
  static void WriteWithXmf(
      const FieldCell<T>& fc, const std::string& fieldname,
      const std::string& rawpath, M& m);

  // Dumps a multidimensional array to `rawpath` and writes metadata
  // to an XMF file in the same directory.
  // The mesh is uniform of size `meshsize` and
  // extends between `xlower` and `xupper`.
  // The domain will have extent 1 and origin at 0.
  template <class T>
  static void WritePlainArrayWithXmf(
      const std::string& rawpath, const std::string& fieldname, const T* data,
      MIdx meshsize, Vect xlower, Vect xupper);

  template <class T>
  static void Read(
      FieldCell<T>& fc, const Meta& meta, const std::string& path, M& m);
};

} // namespace dump
