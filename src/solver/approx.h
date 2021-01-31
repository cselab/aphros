// Created by Petr Karnakov on 15.05.2018
// Copyright 2018 ETH Zurich

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

// sc: scheme:
//   - fou: first order upwind
//   - cd: central differences (mean value)
//   - sou: second order upwind
//   - quick: QUICK
template <class Scal>
std::array<Scal, 3> GetCoeff(ConvSc);

// Linear extrapolation.
// xt: target
// x0,x1: points
// v0,v1: values
template <class T, class Scal>
T UExtrap(Scal xt, Scal x0, const T& v0, Scal x1, const T& v1) {
  return v0 + (v1 - v0) * ((xt - x0) / (x1 - x0));
}

// Smoothens fieldcell with node-based averaging.
// fc: fieldcell [s]
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void SmoothenNode(FieldCell<T>& fc, M& m, size_t rep);

// Smoothens node field.
// fc: fieldcell [s]
// iters: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void SmoothenNode(FieldNode<T>& fc, M& m, size_t iters);

// Coefficients for approximation of gradient with polynomial.
// x: target point
// z: stencil points
// Output:
// k: such that grad(x) = sum_i (ki * f(zi))
template <class Scal>
std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z);

// Returns GetGradCoeffs(x,z[b:]) preceeded by b zeros.
template <class Scal>
std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z, size_t b);

// Apply boundary conditions to halo cells
template <class T, class M>
void BcApply(FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);

// Apply reflection on all boundaries
// fill: value for other types that CondFaceReflect
template <class T, class M>
void BcReflectAll(FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);

template <class Scal>
Scal Superbee(Scal p, Scal q) {
  if (p > 0. && q > 0.) {
    return std::max(std::min(2 * p, q), std::min(p, 2 * q));
  } else if (p < 0. && q < 0.) {
    return -std::max(std::min(-2 * p, -q), std::min(-p, -2 * q));
  }
  return 0.;
}
