// Created by Petr Karnakov on 30.11.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "cond.h"
#include "embed.h"
#include "geom/mesh.h"
#include "solver.h"

template <class MEB>
struct Approx2 {
  using M = typename MEB::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  template <class T>
  using FieldFaceb = typename EmbedTraits<MEB>::template FieldFaceb<T>;

  // Returns divergence of face field
  // defined as sum of fluxes divided by volume of regular cell.
  static FieldCell<Scal> GetRegularDivergence(
      const FieldFaceb<Scal>& fev, const MEB& eb) {
    auto& m = eb.GetMesh();
    FieldCell<Scal> fcdiv(m, 0);
    for (auto c : eb.Cells()) {
      Scal div = 0;
      eb.LoopNci(c, [&](auto q) {
        const auto cf = eb.GetFace(c, q);
        div += fev[cf] * eb.GetOutwardFactor(c, q);
      });
      fcdiv[c] = div / m.GetVolume(c);
    }
    return fcdiv;
  }
};
