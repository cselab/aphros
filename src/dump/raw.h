// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <string>
#include <ostream>
#include <istream>
#include "geom/mesh.h"

namespace dump {

template <class M>
class Raw {
 public:
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

  enum class Type { UInt16, Float32, Float64 };
  static std::string TypeToString(Type);
  static int GetPrecision(Type);
  static Type StringToType(std::string);
  struct Meta {
    std::string name{"u"};
    std::string binpath;
    Type type{Type::Float64};
    MIdx dimensions; // full size, number of cells
    MIdx start{0}; // hyperslab start
    MIdx stride{1}; // hyperslab stride
    MIdx count; // hypreslab size
    IntIdx seek{0}; // number of bytes to skip
    Vect origin{0}; // position of the grid node with zero index
    Vect spacing{1}; // spacing between grid nodes
  };

  template <class T>
  static void Write(
      const FieldCell<T>& fc, const Meta& meta, std::string path, M& m);
  template <class T>
  static void Read(FieldCell<T>& fc, const Meta& meta, std::string path, M& m);

  static std::string GetXmfTemplate();

  static Meta GetMeta(MIdx start, MIdx stride, const M&);
  static void WriteXmf(std::ostream&, const Meta&);
  static void WriteXmf(const std::string& xmfpath, const Meta&);

  static Meta ReadXmf(std::istream&);
  static Meta ReadXmf(const std::string& xmfpath);
};

} // namespace dump
