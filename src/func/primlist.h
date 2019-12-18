#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "geom/vect.h"

template <class Scal>
struct GPrimitive {
  static constexpr size_t dim = 3;
  using Vect = GVect<Scal, dim>;
  Vect c; // center
  Vect r; // axes in coordinate directions
  Vect n; // normal
  Scal th; // ring thickness
  Scal magn; // magnitude
  std::function<Scal(const GVect<Scal, 3>&)> ls; // level-set
  std::function<Scal(const Rect<GVect<Scal, 3>>&)> inter; // true if intersects
  std::function<Vect(const GVect<Scal, 3>&)> vel; // velocity
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
    // remove comments
    sv = sv.substr(0, sv.find('#', 0));
    // check empty string
    if (std::all_of(
            sv.cbegin(), sv.cend(), [](char c) { return std::isspace(c); })) {
      return std::map<std::string, Scal>();
    }
    std::string pre; // name
    std::stringstream vv(sv);
    vv >> std::skipws;
    {
      std::stringstream tt(sv);
      Scal a;
      tt >> std::skipws;
      tt >> a;
      if (tt.good()) { // first is number
        pre = "";
      } else {
        vv >> pre;
      }
    }
    if (name != pre) {
      return std::map<std::string, Scal>();
    }
    std::stringstream kk(sk);
    std::map<std::string, Scal> d;
    std::string k;
    Scal v;
    while (true) {
      vv >> v;
      kk >> k;
      if (vv && kk) {
        d[k] = v;
      } else {
        break;
      }
    }
    if (d.size() < n && kk) {
      throw std::runtime_error(
          "PrimList: missing field '" + k + "' in '" + sv +
          "' while parsing keys '" + sk + "'");
    }
    if (d.size() < n && !kk) {
      throw std::runtime_error(
          "PrimList: no keys for " + std::to_string(n) +
          " required fields in '" + sk + "'");
    }
    if (vv && !kk) {
      throw std::runtime_error(
          "PrimList: no key for '" + std::to_string(v) + "' in '" + sv +
          "' while parsing keys '" + sk + "'");
    }
    return d;
  }

  static std::vector<Primitive> Parse(std::istream& f, bool verb, size_t edim) {
    std::vector<Primitive> pp;

    f >> std::skipws;

    // default
    auto def = [](std::map<std::string, Scal>& d, std::string k, Scal a) {
      if (!d.count(k)) {
        d[k] = a;
      }
    };

    // Read until eof
    while (f) {
      std::string s;
      std::getline(f, s);
      std::map<std::string, Scal> d;

      // sphere
      d = Parse("s", s, "x y z rx ry rz", 4);
      if (d.empty()) {
        d = Parse("", s, "x y z rx ry rz", 4);
      }

      if (!d.empty()) {
        def(d, "ry", d["rx"]);
        def(d, "rz", d["ry"]);

        Primitive p;
        p.c[0] = d["x"];
        p.c[1] = d["y"];
        p.c[2] = d["z"];
        p.r[0] = d["rx"];
        p.r[1] = d["ry"];
        p.r[2] = d["rz"];

        p.ls = [edim, p](const Vect& x) -> Scal {
          Vect xd = (x - p.c) / p.r;
          if (edim == 2) {
            xd[2] = 0.;
          }
          return (1. - xd.sqrnorm()) * sqr(p.r.min());
        };

        p.inter = [edim, p](const Rect<Vect>& rect) -> bool {
          const Rect<Vect> rectbig(rect.lb - p.r, rect.rt + p.r);
          return rectbig.IsInside(p.c);
        };

        pp.push_back(p);
      }

      // ring
      d = Parse("ring", s, "x y z nx ny nz r th", 8);
      if (!d.empty()) {
        Primitive p;
        p.c[0] = d["x"];
        p.c[1] = d["y"];
        p.c[2] = d["z"];
        p.r[0] = d["r"];
        p.n[0] = d["nx"];
        p.n[1] = d["ny"];
        p.n[2] = d["nz"];
        p.th = d["th"];
        p.n /= p.n.norm();

        p.ls = [p](const Vect& x) -> Scal {
          Vect d = x - p.c;
          Scal xn = d.dot(p.n);
          Scal xt = (d - p.n * xn).norm();
          Scal r = p.r[0];
          return sqr(p.th) - (sqr(xn) + sqr(xt - r));
        };

        p.inter = [edim, p](const Rect<Vect>&) -> bool { return true; };

        pp.push_back(p);
      }

      // box
      d = Parse("box", s, "x y z rx ry rz", 4);
      if (!d.empty()) {
        def(d, "ry", d["rx"]);
        def(d, "rz", d["ry"]);

        Primitive p;
        p.c[0] = d["x"];
        p.c[1] = d["y"];
        p.c[2] = d["z"];
        p.r[0] = d["rx"];
        p.r[1] = d["ry"];
        p.r[2] = d["rz"];

        p.ls = [edim, p](const Vect& x) -> Scal {
          Vect xd = (x - p.c) / p.r;
          if (edim == 2) {
            xd[2] = 0.;
          }
          return (1. - xd.max());
        };

        p.inter = [edim, p](const Rect<Vect>&) -> bool { return true; };

        pp.push_back(p);
      }
    }

    if (verb) {
      std::cout << "Read " << pp.size() << " primitives." << std::endl;
    }
    return pp;
  }
  static std::vector<Primitive> ParseVel(
      std::string fn, bool verb, size_t edim) {
    (void)dim;
    (void)edim;
    std::vector<Primitive> pp;
    std::ifstream f(fn);
    if (!f.good() && verb) {
      throw std::runtime_error("Can't open primitive list '" + fn + "'");
    }

    f >> std::skipws;

    // default
    auto def = [](std::map<std::string, Scal>& d, std::string k, Scal a) {
      if (!d.count(k)) {
        d[k] = a;
      }
    };
    (void)def;

    // Read until eof
    while (f) {
      std::string s;
      std::getline(f, s);
      std::map<std::string, Scal> d;

      // ring
      d = Parse("ring", s, "x y z nx ny nz r th magn", 9);
      if (!d.empty()) {
        Primitive p;
        p.c[0] = d["x"];
        p.c[1] = d["y"];
        p.c[2] = d["z"];
        p.r[0] = d["r"];
        p.n[0] = d["nx"];
        p.n[1] = d["ny"];
        p.n[2] = d["nz"];
        p.th = d["th"];
        p.magn = d["magn"];
        p.n /= p.n.norm();

        p.vel = [p](const Vect& x) -> Vect {
          const Scal eps = 1e-10;
          Vect d = x - p.c;
          Scal xn = d.dot(p.n);
          Scal xt = (d - p.n * xn).norm();
          Scal s2 = sqr(xn) + sqr(xt - p.r[0]);
          Scal sig2 = sqr(p.th);
          // unit radial along plane
          Vect et = (d - p.n * xn) / std::max(eps, xt);
          // unit along circle
          Vect es = p.n.cross(et);
          Scal om = p.magn / (M_PI * sig2) * std::exp(-s2 / sig2);
          return es * om;
        };

        pp.push_back(p);
      }
    }

    if (verb) {
      std::cout << "Read " << pp.size() << " primitives from "
                << "'" << fn << "'" << std::endl;
    }
    return pp;
  }
};

template <class Scal>
std::ostream& operator<<(std::ostream& o, const GPrimitive<Scal>& p) {
  o << "c=" << p.c << " r=" << p.r << " n=" << p.n << " th=" << p.r
    << " inv=" << p.inv;
  return o;
}
