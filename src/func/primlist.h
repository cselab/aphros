// Created by Petr Karnakov on 22.09.2019
// Copyright 2019 ETH Zurich

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
  std::function<Scal(const Vect&)> ls; // level-set
  std::function<bool(const Rect<Vect>&)> inter; // true if intersects rectangle
  std::function<Vect(const Vect&)> vel; // velocity

  Vect c; // center XXX adhoc for GetSphereOverlap
  Vect r; // radius XXX adhoc for GetSphereOverlap

  GPrimitive() {
    inter = [](const Rect<Vect>&) -> bool { return true; };
  }
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

    // sets default
    auto setdef = [](std::map<std::string, Scal>& d, std::string k, Scal a) {
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
        d = Parse("sphere", s, "x y z rx ry rz", 4);
      }
      if (d.empty()) {
        d = Parse("", s, "x y z rx ry rz", 4);
      }

      if (!d.empty()) {
        setdef(d, "ry", d["rx"]);
        setdef(d, "rz", d["ry"]);

        Primitive p;
        const Vect xc(d["x"], d["y"], d["z"]); // center
        const Vect r(d["rx"], d["ry"], d["rz"]); // radius (semi-axes)
        p.c = xc;
        p.r = r;

        p.ls = [edim, xc, r](const Vect& x) -> Scal {
          Vect xd = (x - xc) / r;
          if (edim == 2) {
            xd[2] = 0.;
          }
          return (1. - xd.sqrnorm()) * sqr(r.min());
        };

        p.inter = [xc, r](const Rect<Vect>& rect) -> bool {
          const Rect<Vect> rectbig(rect.lb - r, rect.rt + r);
          return rectbig.IsInside(xc);
        };

        pp.push_back(p);
      }

      // ring
      d = Parse("ring", s, "x y z nx ny nz r th", 8);
      if (!d.empty()) {
        Primitive p;

        const Vect xc(d["x"], d["y"], d["z"]); // center
        Vect n(d["nx"], d["ny"], d["nz"]); // normal
        n /= n.norm();
        const Scal r = d["r"]; // radius
        const Scal th = d["th"]; // thickness

        p.ls = [xc, n, r, th](const Vect& x) -> Scal {
          const Vect dx = x - xc;
          const Scal xn = dx.dot(n);
          const Scal xt = (dx - n * xn).norm();
          return sqr(th) - (sqr(xn) + sqr(xt - r));
        };

        pp.push_back(p);
      }

      // box
      d = Parse("box", s, "x y z rx ry rz", 4);
      if (!d.empty()) {
        setdef(d, "ry", d["rx"]);
        setdef(d, "rz", d["ry"]);

        Primitive p;
        const Vect xc(d["x"], d["y"], d["z"]); // center
        const Vect r(d["rx"], d["ry"], d["rz"]); // radius (semi-axes)

        p.ls = [edim, xc, r](const Vect& x) -> Scal {
          Vect dx = (x - xc) / r;
          if (edim == 2) {
            dx[2] = 0.;
          }
          return 1. - dx.abs().max();
        };

        p.inter = [p](const Rect<Vect>&) -> bool { return true; };

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
    auto setdef = [](std::map<std::string, Scal>& d, std::string k, Scal a) {
      if (!d.count(k)) {
        d[k] = a;
      }
    };
    (void)setdef;

    // Read until eof
    while (f) {
      std::string s;
      std::getline(f, s);
      std::map<std::string, Scal> d;

      // ring
      d = Parse("ring", s, "x y z nx ny nz r th magn", 9);
      if (!d.empty()) {
        Primitive p;
        const Vect xc(d["x"], d["y"], d["z"]); // center
        Vect n(d["nx"], d["ny"], d["nz"]); // normal
        n /= n.norm();
        const Scal r = d["r"]; // radius
        const Scal th = d["th"]; // thickness
        const Scal magn = d["magn"]; // magnitude

        p.vel = [xc, n, r, th, magn](const Vect& x) -> Vect {
          const Scal eps = 1e-10;
          const Vect dx = x - xc;
          const Scal xn = dx.dot(n);
          const Scal xt = (dx - n * xn).norm();
          const Scal s2 = sqr(xn) + sqr(xt - r);
          const Scal sig2 = sqr(th);
          // unit radial along plane
          const Vect et = (dx - n * xn) / std::max(eps, xt);
          // unit along circle
          const Vect es = n.cross(et);
          const Scal om = magn / (M_PI * sig2) * std::exp(-s2 / sig2);
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
