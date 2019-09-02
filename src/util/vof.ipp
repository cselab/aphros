#pragma once

#include <vector>
#include <memory>
#include <limits>

#include "vof.h"
#include "geom/mesh.h"
#include "dump/vtk.h"
#include "solver/reconst.h"

template <class M_>
struct UVof<M_>::Imp {
  using M = M_;
  using R = Reconst<Scal>;
  static constexpr size_t dim = M::dim;

  Imp(M& m) : m(m) {}

  M& m;
  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<Scal> dlc_; // dump poly
};

template <class M_>
UVof<M_>::UVof(M& m) : imp(new Imp(m)) {}

template <class M_>
UVof<M_>::~UVof() = default;

template <class M_>
void UVof<M_>::Dump() {
  return;
}
