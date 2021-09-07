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
  using Meta = typename Xmf<Vect>::Meta;

  template <class T>
  static void Write(
      const std::string& path, const std::vector<MIdx>& starts,
      const std::vector<MIdx>& sizes, const std::vector<std::vector<T>>& data,
      MIdx global_size, Type type, const MpiWrapper& mpi);
  template <class T>
  static void Write(
      const FieldCell<T>& fc, const Meta& meta, std::string path, M& m);
  // Writes field to raw-file and metadata to xmf-file.
  // fc: field to write
  // field_name: dataset name to put in metadata
  // raw_path: output path to file with extension `.raw`
  template <class T>
  static void WriteWithXmf(
      const FieldCell<T>& fc, std::string fieldname, std::string rawpath, M& m);
  template <class T>
  static void Read(FieldCell<T>& fc, const Meta& meta, std::string path, M& m);
};

} // namespace dump
