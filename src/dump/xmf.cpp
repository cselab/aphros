// Created by Petr Karnakov on 09.02.2021
// Copyright 2021 ETH Zurich

#include "geom/mesh.h"
#include "util/logger.h"

#include "xmf.ipp"

namespace dump {

Type StringToType(std::string s) {
  Type res;
  if (s == "UShort") {
    res = Type::UInt16;
  } else if (s == "Float") {
    res = Type::Float32;
  } else if (s == "Double") {
    res = Type::Float64;
  } else {
    fassert(false, "Unknown type=" + s);
  }
  return res;
}

std::string TypeToString(Type type) {
  switch (type) {
    case Type::UInt16:
      return "UShort";
    case Type::Float32:
      return "Float";
    case Type::Float64:
      return "Double";
    default:
      fassert(false);
  }
  return "";
}

int GetPrecision(Type type) {
  switch (type) {
    case Type::UInt16:
      return 2;
    case Type::Float32:
      return 4;
    case Type::Float64:
      return 8;
    default:
      fassert(false);
  }
  return 8;
}

#define X(dim) template class Xmf<generic::Vect<double, dim>>;
X(1)
X(2)
X(3)
X(4)
#undef X

#define XX(M)                                   \
  template typename Xmf<typename M::Vect>::Meta \
  Xmf<typename M::Vect>::GetMeta<M>(            \
      const M&, typename M::MIdx, typename M::MIdx);
#define COMMA ,

#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX

#undef X
#undef XX

} // namespace dump
