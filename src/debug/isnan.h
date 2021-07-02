// Created by Petr Karnakov on 03.10.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <exception>
#include <limits>
#include <sstream>
#include <string>

#include "geom/field.h"
#include "geom/range.h"
#include "geom/vect.h"
#include "solver/embed.h"
#include "util/format.h"
#include "util/logger.h"

template <class T, class Idx>
bool IsFinite(const GField<T, Idx>& u) {
  for (auto i : u.GetRange()) {
    if (!IsFinite(u[i])) {
      return false;
    }
  }
  return true;
}

template <class T, class Idx>
bool IsNan(const GField<T, Idx>& u) {
  for (auto i : u.GetRange()) {
    if (IsNan(u[i])) {
      return true;
    }
  }
  return false;
}

// Checks for Nan in field.
// u: field
// n: name
// r: range
// w: where, name of location (e.g. __FILE__)
template <class T, class Idx>
bool CheckNan(
    const GField<T, Idx>& u, std::string n, const GRange<Idx>& range,
    std::string msg = "") {
  if (msg != "") {
    msg += ": ";
  }
  for (auto i : range) {
    fassert(
        !IsNan(u[i]),
        util::Format("{} Nan field {} at i={}", msg, n, size_t(i)));
  }
  return false;
}

// Checks for Nan in inner cells/faces.
// u: field
// n: name
// m: mesh
// w: where, name of location (e.g. __FILE__)
template <class T, class Idx, class M>
bool CheckNan(
    const GField<T, Idx>& u, std::string n, const M& m, std::string msg = "") {
  if (msg != "") {
    msg += ": ";
  }
  for (auto i : m.template GetRangeIn<Idx>()) {
    fassert(
        !IsNan(u[i]),
        util::Format("{} Nan field {} at x={}", msg, n, m.GetCenter(i)));
  }
  return false;
}

template <class T, class M>
bool CheckNan(
    const FieldEmbed<T>& u, std::string n, const M& m, std::string w = "") {
  return CheckNan(u.GetFieldCell(), n, m, w) ||
         CheckNan(u.GetFieldFace(), n, m, w);
}

// Returns path of file relative to src/
inline std::string SrcPath(const std::string& s) {
  std::string k = "src/";
  size_t i = s.rfind(k);
  if (i != std::string::npos) {
    return s.substr(i + k.length());
  }
  return s;
}

// CheckNan for field if condition evaluates to true.
// field: instance of GField
// cond: bool expression
// Requires in scope:
// m: mesh or GRange
#define CHECKNAN(field, cond)                                         \
  if (cond) {                                                         \
    CheckNan(                                                         \
        field, #field, m,                                             \
        SrcPath(__FILE__) + ":" + std::to_string(__LINE__) + " in " + \
            std::string(__func__));                                   \
  }
