// Created by Petr Karnakov on 09.02.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <iosfwd>
#include <string>

#include "geom/idx.h"

namespace dump {

enum class Type { UInt16, Float32, Float64 };
std::string TypeToString(Type);
Type StringToType(std::string);
int GetPrecision(Type);

template <class Vect_>
class Xmf {
 public:
  using Vect = Vect_;
  using Scal = typename Vect::Scal;
  using MIdx = generic::MIdx<Vect::dim>;
  static constexpr size_t dim = Vect::dim;

  static std::string GetXmfTemplate();

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

  template <class M>
  static Meta GetMeta(const M& m, MIdx start = MIdx(0), MIdx stride = MIdx(1));

  static void WriteXmf(std::ostream&, const Meta&);
  static void WriteXmf(const std::string& xmfpath, const Meta&);

  static Meta ReadXmf(std::istream&);
  static Meta ReadXmf(const std::string& xmfpath);
};

} // namespace dump
