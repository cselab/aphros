// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <string>
#include "geom/mesh.h"

template <class M>
class Hdf {
 public:
  using Scal = typename M::Scal;
  using MIdx = typename M::MIdx;
  static void Write(
      const FieldCell<typename M::Scal>& fc, std::string path, M& m);
};
