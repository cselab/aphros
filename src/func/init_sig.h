// Created by Petr Karnakov on 24.10.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "geom/block.h"
#include "geom/field.h"
#include "geom/vect.h"
#include "parse/vars.h"
#include "solver/reconst.h"

// Volume fraction field.
// par: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc,const M& m)>
// fc: field to fill [i]
// m: mesh
template <class M>
std::function<void(FieldCell<typename M::Scal>&, const M&)> CreateInitSig(
    const Vars& par) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  std::function<void(FieldCell<Scal>&, const M&)> g; // result

  std::string v = par.String["init_sig"];
  if (v == "uniform") {
    Scal sig = par.Double["sigma"];
    g = [sig](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        fc[c] = sig;
      }
    };
  } else if (v == "ge") {
    Scal sig = par.Double["sigma"];
    Scal k = par.Double["sig_k"];
    g = [sig, k](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        fc[c] = sig * std::max(1. - k * std::abs(x[0] - 0.5), 0.1);
      }
    };
  } else if (v == "linear") {
    Scal sig = par.Double["sigma"];
    Vect x0(par.Vect["sig_x"]);
    Vect grad(par.Vect["sig_grad"]);
    g = [sig, x0, grad](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        fc[c] = sig + grad.dot(x - x0);
      }
    };
  } else if (v == "linearlow") {
    Scal sig = par.Double["sigma"];
    Vect x0(par.Vect["sig_x"]);
    Vect grad(par.Vect["sig_grad"]);
    g = [sig, x0, grad](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        fc[c] = sig + std::min(0., grad.dot(x - x0));
      }
    };
  } else if (v == "linearhigh") {
    Scal sig = par.Double["sigma"];
    Vect x0(par.Vect["sig_x"]);
    Vect grad(par.Vect["sig_grad"]);
    g = [sig, x0, grad](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        fc[c] = sig + std::max(0., grad.dot(x - x0));
      }
    };
  } else {
    fassert(false, "Unknown init_sig=" + v);
  }
  return g;
}
