#pragma once

#include <limits>
#include <exception>
#include <sstream>
#include <string>

#include "geom/vect.h"
#include "geom/field.h"
#include "geom/range.h"

template <class Scal>
bool IsFinite(Scal a) {
  return std::isfinite(a);
}

template <class T, class Idx>
bool IsFinite(const GField<T, Idx>& u) {
  for (auto i : u.GetRange()) {
    if (!IsFinite(u[i])) {
      return false;
    }
  }
  return true;
}

template <class Scal>
bool IsNan(Scal a) {
  return std::isnan(a);
}

template <class Scal, size_t dim>
bool IsNan(const GVect<Scal, dim>& a) {
  for (size_t i = 0; i < dim; ++i) {
    if (IsNan(a[i])) {
      return true;
    }
  }
  return false;
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
bool CheckNan(const GField<T, Idx>& u, std::string n, 
              const GRange<Idx>& r, std::string w="") {
  if (w != "") {
    w += ": ";
  }
  for (auto i : r) {
    if (IsNan(u[i])) {
      std::stringstream s;
      s << w << "Nan " << n << " at i=" << size_t(i);
      throw std::runtime_error(s.str());
    }
  }
  return false;
}

// Checks for Nan in inner cells/faces.
// u: field
// n: name
// m: mesh
// w: where, name of location (e.g. __FILE__)
template <class T, class Idx, class M>
bool CheckNan(const GField<T, Idx>& u, std::string n, 
              const M& m, std::string w="") {
  if (w != "") {
    w += ": ";
  }
  for (auto i : m.template GetIn<Idx>()) {
    if (IsNan(u[i])) {
      std::stringstream s;
      s << w << "Nan " << n << " at x=" << m.GetCenter(i);
      throw std::runtime_error(s.str());
    }
  }
  return false;
}

// CheckNan for field with file and line.
// F: instance of GField
// Requires in scope:
// m: mesh or GRange
#define CHECKNAN(F) \
  CheckNan(F, #F, m, std::string(__FILE__) + ":" + std::to_string(__LINE__));
