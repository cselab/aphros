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
  template <class T>
  static void Read(FieldCell<T>& fc, const Meta& meta, std::string path, M& m);
};

} // namespace dump
