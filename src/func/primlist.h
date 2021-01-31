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

template <class Vect_>
struct Primitive {
  using Vect = Vect_;
  using Scal = typename Vect::Scal;
  static constexpr size_t dim = Vect::dim;

  Primitive(std::string name_)
      : name(name_), inter([](const Rect<Vect>&) -> bool { return true; }) {}
  Primitive() : Primitive("") {}
  friend std::ostream& operator<<(std::ostream& o, const Primitive<Vect>& p) {
    o << "name='" << p.name << "'";
    o << " mod='" << (p.mod_minus ? "-" : "") << (p.mod_and ? "&" : "") << "'";
    return o;
  }

  std::string name;
  bool mod_minus = false;
  bool mod_and = false;
  std::function<Scal(const Vect&)> ls; // level-set
  std::function<bool(const Rect<Vect>&)> inter; // true if intersects rectangle
  std::function<Vect(const Vect&)> velocity; // velocity

  Vect c; // center XXX adhoc for GetSphereOverlap
  Vect r; // radius XXX adhoc for GetSphereOverlap
};

template <class Vect_>
struct VelocityPrimitive {
  using Vect = Vect_;
  using Scal = typename Vect::Scal;
  static constexpr size_t dim = Vect::dim;

  VelocityPrimitive() : inter([](const Rect<Vect>&) { return true; }) {}
  friend std::ostream& operator<<(
      std::ostream& o, const VelocityPrimitive<Vect>& p) {
    o << "name='" << p.name << "'";
    o << " mod='" << (p.mod_minus ? "-" : "") << (p.mod_and ? "&" : "") << "'";
    return o;
  }

  std::function<Vect(const Vect&)> velocity; // velocity
  // true if rectangle contains points with non-zero velocity
  std::function<bool(const Rect<Vect>&)> inter;
};

} // namespace generic

template <class Vect_>
struct UPrimList {
  using Vect = Vect_;
  using Scal = typename Vect::Scal;
  static constexpr size_t dim = Vect::dim;

  using Primitive = generic::Primitive<Vect>;
  using VelocityPrimitive = generic::VelocityPrimitive<Vect>;

  // Parses a list of primitives in stream buf.
  // edim: effective dimension, 2 or 3 (ignores z-component if edim=2)
  static std::vector<Primitive> GetPrimitives(std::istream& buf, size_t edim);
  static std::vector<VelocityPrimitive> GetVelocityPrimitives(
      std::istream& buf, size_t edim);
};
