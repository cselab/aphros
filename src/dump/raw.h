// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <string>
#include "geom/mesh.h"

namespace dump {

template <class M>
class Raw {
 public:
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

  enum class Type { UInt16, Float32, Float64 };
  struct Meta {
    IntIdx seek{0};
    MIdx size;
    // hyperslab
    MIdx start{0};
    MIdx stride{1};
    MIdx count;
    Type type;
  };

  template <class T>
  static void Write(
      const FieldCell<T>& fc, const Meta& meta, std::string path, M& m);
  template <class T>
  static void Read(FieldCell<T>& fc, const Meta& meta, std::string path, M& m);

  static void WriteXmf(
      std::string xmfpath, std::string name, Type type, Vect origin,
      Vect spacing, MIdx dims, std::string rawpath);
  static void WriteXmf(
      std::string xmfpath, std::string name, Type type, std::string rawpath,
      const M& m);
};

} // namespace dump
