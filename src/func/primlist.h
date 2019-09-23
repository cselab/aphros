#pragma once

#include <cmath>
#include <stdexcept>
#include <functional>
#include <limits>
#include <sstream>
#include <ostream>
#include <vector>
#include <map>
#include <fstream>

#include "geom/vect.h"

template <class Scal>
struct GPrimitive {
  enum class Type { sphere, ring };
  static constexpr size_t dim = 3;
  using Vect = GVect<Scal, dim>;
  Type type;
  Vect c;  // center
  Vect r;  // axes in coordinate directions
  Vect n;  // normal
  Scal th; // ring thickness
  bool inv; // flip
  std::function<Scal(const GVect<Scal, 3>&)> ls; // level-set
};


template <class Scal>
struct UPrimList {
  static constexpr size_t dim = 3;
  using Vect = GVect<Scal, dim>;
  using Primitive = GPrimitive<Scal>;

  // Parses string of format:
  // <name> <r[k1]> <r[k2]> ...
  // and returns map r[k] if name matches, otherwise returns empty map.
  // name: prefix
  // sv: list of values
  // sk: list of keys
  // n: number of required fields
  static std::map<std::string, Scal> Parse(
      std::string name, std::string sv, std::string sk, size_t n) {
    // remove comment
    sv = sv.substr(0, sv.find('#', 0));
    // check empty string
    if (std::all_of(sv.cbegin(), sv.cend(),
                    [](char c){ return std::isspace(c); })) {
      return std::map<std::string, Scal>();
    }
    std::stringstream vv(sv);
    vv >> std::skipws;
    if (name != "") {
      std::string q;
      vv >> q;
      if (q != name) {
        return std::map<std::string, Scal>();
      }
    }
    std::stringstream kk(sk);
    std::map<std::string, Scal> r;
    std::string k;
    Scal v;
    while (true) {
      vv >> v;
      kk >> k;
      if (vv && kk) {
        r[k] = v;
      } else {
        break;
      }
    }
    if (r.size() < n && kk) {
      throw std::runtime_error(
          "PrimList: missing field '" + k + "' in '" + sv +
          "' while parsing keys '" + sk + "'");
    }
    if (r.size() < n && !kk) {
      throw std::runtime_error(
          "PrimList: no keys for " +std::to_string(n) + 
          " required fields in '" + sk + "'");
    }
    if (vv && !kk) {
      throw std::runtime_error(
          "PrimList: no key for '" + std::to_string(v) + "' in '" + sv +
          "' while parsing keys '" + sk + "'");
    }
    return r;
  }

  static std::vector<Primitive> Parse(
      std::string fn, bool verb, size_t edim) {
    std::vector<Primitive> pp;
    std::ifstream f(fn);
    if (!f.good() && verb) {
      throw std::runtime_error("Can't open particle list '" + fn + "'");
    }

    f >> std::skipws;

    // default
    auto def = [](std::map<std::string, Scal>& r, std::string k, Scal d) {
      if (!r.count(k)) { r[k] = d; }
    };

    // Read until eof
    while (f) {
      std::string s;
      std::getline(f, s);
      std::map<std::string, Scal> r;

      r  = Parse("s", s, "x y z rx ry rz", 4);
      if (r.empty()) {
        r = Parse("", s, "x y z rx ry rz", 4);
      }

      if (!r.empty()) { // sphere
        def(r, "ry", r["rx"]);
        def(r, "rz", r["ry"]);

        Primitive p;
        p.c[0] = r["x"];
        p.c[1] = r["y"];
        p.c[2] = r["z"];
        p.r[0] = r["rx"];
        p.r[1] = r["ry"];
        p.r[2] = r["rz"];

        using Type = typename Primitive::Type;
        p.type = Type::sphere;

        p.ls = [edim,p](const Vect& x) -> Scal {
          Vect xd = (x - p.c) / p.r;
          if (edim == 2) {
            xd[2] = 0.;
          }
          return 1. - xd.sqrnorm();
        };

        pp.push_back(p);
      }
    }

    if (verb) {
      std::cout << "Read " << pp.size() << " particles from " 
          << "'" << fn << "'" << std::endl;
    }
    return pp;
  }
};

template <class Scal>
std::ostream& operator<<(std::ostream& o,
                         const GPrimitive<Scal>& p) {
  o << "c=" << p.c
    << " r=" << p.r
    << " n=" << p.n
    << " th=" << p.r
    << " inv=" << p.inv;
  return o;
}

