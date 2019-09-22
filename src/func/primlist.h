#pragma once

#include <cmath>
#include <stdexcept>
#include <functional>
#include <limits>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>

#include "geom/vect.h"

template <class Scal>
struct UPrimList {
  static constexpr size_t dim = 3;
  using Vect = GVect<Scal, dim>;

  struct Primitive {
    Vect c;  // center
    Vect r;  // axes in coordinate directions
  };

  // Parses string of format:
  // <name> <r[k1]> <r[k2]> ...
  // and returns map r[k] if name matches, otherwise returns empty map.
  // name: prefix
  // sv: list of values
  // sk: list of keys
  // n: number of required fields
  static std::map<std::string, Scal> Parse(
      std::string name, std::string sv, std::string sk, size_t n) {
    std::map<std::string, Scal> r;
    std::stringstream vv(sv);
    std::stringstream kk(sk);
    if (name != "") {
      std::string q;
      vv >> q;
      if (q != name) {
        return std::map<std::string, Scal>();
      }
    }
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

  static std::vector<Primitive> Parse(std::string fn, bool verb) {
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
      Primitive p;
      std::string s;
      std::getline(f, s);

      std::map<std::string, Scal> r = Parse("s", s, "x y z rx ry rz", 4);

      def(r, "ry", r["rx"]);
      def(r, "rz", r["ry"]);

      p.c[0] = r["x"];
      p.c[1] = r["y"];
      p.c[2] = r["z"];
      p.r[0] = r["rx"];
      p.r[1] = r["ry"];
      p.r[2] = r["rz"];

      pp.push_back(p);
    }

    if (verb) {
      std::cout << "Read " << pp.size() << " particles from " 
          << "'" << fn << "'" << std::endl;
    }
    return pp;
  }
};
