// Created by Petr Karnakov on 22.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <algorithm>
#include <functional>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

#include "geom/vect.h"

namespace generic {

template <class Scal>
struct Primitive {
  Primitive(std::string name_)
      : name(name_), inter([](const Rect<Vect>&) -> bool { return true; }) {}
  Primitive() : Primitive("") {}
  friend std::ostream& operator<<(std::ostream& o, const Primitive<Scal>& p) {
    o << "name='" << p.name << "'";
    o << " mod='" << (p.mod_minus ? "-" : "") << (p.mod_and ? "&" : "") << "'";
    return o;
  }

  static constexpr size_t dim = 3;
  using Vect = generic::Vect<Scal, dim>;

  std::string name;
  bool mod_minus = false;
  bool mod_and = false;
  std::function<Scal(const Vect&)> ls; // level-set
  std::function<bool(const Rect<Vect>&)> inter; // true if intersects rectangle
  std::function<Vect(const Vect&)> velocity; // velocity

  Vect c; // center XXX adhoc for GetSphereOverlap
  Vect r; // radius XXX adhoc for GetSphereOverlap
};

template <class Scal>
struct VelocityPrimitive {
  VelocityPrimitive() : inter([](const Rect<Vect>&) { return true; }) {}
  friend std::ostream& operator<<(
      std::ostream& o, const VelocityPrimitive<Scal>& p) {
    o << "name='" << p.name << "'";
    o << " mod='" << (p.mod_minus ? "-" : "") << (p.mod_and ? "&" : "") << "'";
    return o;
  }

  static constexpr size_t dim = 3;
  using Vect = generic::Vect<Scal, dim>;

  std::function<Vect(const Vect&)> velocity; // velocity
  // true if rectangle contains points with non-zero velocity
  std::function<bool(const Rect<Vect>&)> inter;
};

} // namespace generic

template <class Scal>
struct UPrimList {
  static constexpr size_t dim = 3;
  using Vect = generic::Vect<Scal, dim>;
  using Primitive = generic::Primitive<Scal>;
  using VelocityPrimitive = generic::VelocityPrimitive<Scal>;

  // Parses a list of primitives in stream buf.
  // edim: effective dimension, 2 or 3 (ignores z-component if edim=2)
  static std::vector<Primitive> GetPrimitives(std::istream& buf, size_t edim);
  static std::vector<VelocityPrimitive> GetVelocityPrimitives(
      std::istream& buf, size_t edim);
};
