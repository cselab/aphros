#include <stdexcept>
#include <string>

#include "solver.h"

namespace solver {

std::string GetName(Layers l) {
  switch (l) {
    case Layers::time_curr: { return "time_curr"; }
    case Layers::time_prev: { return "time_prev"; }
    case Layers::iter_curr: { return "iter_curr"; }
    case Layers::iter_prev: { return "iter_prev"; }
    default: {
      throw std::runtime_error("GetName: Unknown layer");
    }
  }
}

ConvSc GetConvSc(std::string s) {
  auto l = {ConvSc::fou, ConvSc::cd, ConvSc::sou, ConvSc::quick};
  for (ConvSc sc : l) {
    if (GetName(sc) == s) {
      return sc;
    }
  }
  throw std::runtime_error("ConvSc: invalid name=" + s);
}

std::string GetName(ConvSc sc) {
  switch (sc) {
    case ConvSc::fou: return "fou";
    case ConvSc::cd: return "cd";
    case ConvSc::sou: return "sou";
    case ConvSc::quick: return "quick";
    default: throw std::runtime_error("GetName: invalid ConvSc");
  }
  return "";
}


} // namespace solver
