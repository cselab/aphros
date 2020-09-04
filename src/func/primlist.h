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

namespace generic {

template <class Scal>
struct Primitive {
  static constexpr size_t dim = 3;
  using Vect = generic::Vect<Scal, dim>;

  std::string name;
  bool mod_minus = false;
  bool mod_and = false;
  std::function<Scal(const Vect&)> ls; // level-set
  std::function<bool(const Rect<Vect>&)> inter; // true if intersects rectangle
  std::function<Vect(const Vect&)> vel; // velocity

  Vect c; // center XXX adhoc for GetSphereOverlap
  Vect r; // radius XXX adhoc for GetSphereOverlap

  Primitive(std::string name)
      : name(name), inter([](const Rect<Vect>&) -> bool { return true; }) {}
  Primitive() : Primitive("") {}

  friend std::ostream& operator<<(std::ostream& o, const Primitive<Scal>& p) {
    o << "name='" << p.name << "'";
    o << " mod='" << (p.mod_minus ? "-" : "") << (p.mod_and ? "&" : "") << "'";
    return o;
  }
};

} // namespace generic

template <class Scal>
struct UPrimList {
  static constexpr size_t dim = 3;
  using Vect = generic::Vect<Scal, dim>;
  using Primitive = generic::Primitive<Scal>;

  static bool ParseSphere(std::string s, size_t edim, Primitive& p) {
    auto d = GetMap("sphere", s, "cx cy cz rx ry rz", 4);
    if (d.empty()) {
      d = GetMap("s", s, "cx cy cz rx ry rz", 4);
    }
    if (d.empty()) {
      return false;
    }
    SetDefault(d, "ry", d["rx"]);
    SetDefault(d, "rz", d["ry"]);

    const Vect xc(d["cx"], d["cy"], d["cz"]); // center
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
    p.name = "sphere";
    return true;
  }
  static bool ParseRing(std::string s, Primitive& p) {
    auto d = GetMap("ring", s, "cx cy cz nx ny nz r th", 8);
    if (d.empty()) {
      return false;
    }
    const Vect xc(d["cx"], d["cy"], d["cz"]); // center
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
    p.name = "ring";
    return true;
  }
  static bool ParseBox(std::string s, size_t edim, Primitive& p) {
    auto d = GetMap("box", s, "cx cy cz rx ry rz rotz", 4);
    if (d.empty()) {
      return false;
    }
    SetDefault(d, "ry", d["rx"]);
    SetDefault(d, "rz", d["ry"]);
    SetDefault(d, "rotz", 0);

    const Vect xc(d["cx"], d["cy"], d["cz"]); // center
    const Vect r(d["rx"], d["ry"], d["rz"]); // radius (semi-axes)
    const Scal rotz = d["rotz"]; // rotation angle around z in degrees
    const Scal rotz_cos = std::cos(-M_PI * rotz / 180.);
    const Scal rotz_sin = std::sin(-M_PI * rotz / 180.);

    p.ls = [edim, xc, r, rotz_cos, rotz_sin](const Vect& x) -> Scal {
      Vect dx = x - xc;
      dx = Vect(
          dx[0] * rotz_cos - dx[1] * rotz_sin,
          dx[0] * rotz_sin + dx[1] * rotz_cos, dx[2]);
      if (edim == 2) {
        dx[2] = 0;
      }
      return (r - dx.abs()).min();
    };
    p.name = "box";
    return true;
  }
  static bool ParseRoundBox(std::string s, size_t edim, Primitive& p) {
    auto d = GetMap("roundbox", s, "cx cy cz rx ry rz round", 7);
    if (d.empty()) {
      return false;
    }

    const Vect xc(d["cx"], d["cy"], d["cz"]); // center
    const Vect r(d["rx"], d["ry"], d["rz"]); // radius (semi-axes)
    const Scal round = d["round"]; // radius of rounding

    p.ls = [edim, xc, r, round](const Vect& x) -> Scal {
      Vect dx = (x - xc).abs();
      if (edim == 2) {
        dx[2] = 0;
      }
      const Vect q = (dx - r + Vect(round)).max(Vect(0));
      return round - q.norm();
    };
    p.name = "roundbox";
    return true;
  }
  static bool ParseSmoothStep(std::string s, Primitive& p) {
    // smooth step (half of diverging channel from almgren1997)
    //  +   +   +   +   +   +  +
    //              -------------       ____t
    //  +   +   +  /c   -   -  -       |
    // ------------                    |n
    //  -   -   -   -   -   -  -
    auto d = GetMap("smooth_step", s, "cx cy cz nx ny nz tx ty tz ln lt", 4);
    if (d.empty()) {
      return false;
    }

    const Vect xc(d["cx"], d["cy"], d["cz"]); // center
    Vect n(d["nx"], d["ny"], d["nz"]); // normal
    n /= n.norm();
    Vect t(d["tx"], d["ty"], d["tz"]); // tangent
    t -= n * n.dot(t);
    t /= t.norm();
    const Scal ln = d["ln"]; // half-length in normal direction
    const Scal lt = d["lt"]; // half-length in tangential direction

    p.ls = [xc, n, t, ln, lt](const Vect& x) -> Scal {
      const Vect dx = x - xc;
      const Scal dn = n.dot(dx) / ln;
      const Scal dt = t.dot(dx) / lt;
      const Scal q = dt < -1 ? 1 : dt > 1 ? -1 : -std::sin(M_PI * 0.5 * dt);
      return (q - dn) * lt;
    };
    p.name = "smooth_step";
    return true;
  }
  static bool ParseCylinder(std::string s, size_t edim, Primitive& p) {
    // c: center
    // t: axis
    // r: radius
    // t0,t1: length between center
    auto d = GetMap("cylinder", s, "cx cy cz tx ty tz r t0 t1 ", 4);
    if (d.empty()) {
      return false;
    }

    const Vect xc(d["cx"], d["cy"], d["cz"]);
    Vect t(d["tx"], d["ty"], d["tz"]);
    t /= t.norm();
    const Scal r = d["r"];
    const Scal t0 = d["t0"];
    const Scal t1 = d["t1"];

    p.ls = [xc, t, r, t0, t1, edim](const Vect& x) -> Scal {
      Vect dx = x - xc;
      if (edim == 2) {
        dx[2] = 0;
      }
      const Scal dt = t.dot(dx);
      const Scal dr = (dx - t * dt).norm();
      Scal q = r - dr;
      if (dt < t0) q = std::min(q, dt - t0);
      if (dt > t1) q = std::min(q, t1 - dt);
      return q;
    };
    p.name = "cylinder";
    return true;
  }

  // Parses a list of primitives in stream f.
  // verb: verbosity (debug messages to std::cout)
  // edim: effective dimension, 2 or 3 (ignores z-component if edim=2)
  static std::vector<Primitive> Parse(std::istream& f, bool verb, size_t edim) {
    std::vector<Primitive> pp;

    f >> std::skipws;

    // Read until eof
    while (f) {
      std::string s;
      std::getline(f, s);
      s = RemoveComment(s);
      if (IsEmpty(s)) {
        continue;
      }

      if (IsNumber(GetWord(s, 0))) {
        s = "sphere " + s;
      }

      Primitive p;

      while (true) {
        auto pair = ExtractModifier(s, "-&");
        if (pair.first == ' ') {
          break;
        }
        s = pair.second;
        switch (pair.first) {
          case '-':
            p.mod_minus = true;
            break;
          case '&':
            p.mod_and = true;
            break;
          default:
            throw std::runtime_error(
                std::string("PrimList: unknown mod='") + pair.first + "'");
        }
      }

      bool r = false;
      if (!r) r = ParseSphere(s, edim, p);
      if (!r) r = ParseRing(s, p);
      if (!r) r = ParseBox(s, edim, p);
      if (!r) r = ParseRoundBox(s, edim, p);
      if (!r) r = ParseSmoothStep(s, p);
      if (!r) r = ParseCylinder(s, edim, p);

      if (p.mod_minus) {
        p.inter = [](const Rect<Vect>&) -> bool { return true; };
      }

      if (!r) {
        throw std::runtime_error(
            "PrimList: primitive not recognized in '" + s + "'");
      }
      pp.push_back(p);
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

    // Read until eof
    while (f) {
      std::string s;
      std::getline(f, s);
      std::map<std::string, Scal> d;

      // ring
      d = GetMap("ring", s, "cx cy cz nx ny nz r th magn", 9);
      if (!d.empty()) {
        Primitive p;
        const Vect xc(d["cx"], d["cy"], d["cz"]); // center
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
      // ring
      d = GetMap("gauss2d", s, "cx cy sig magn", 4);
      if (!d.empty()) {
        Primitive p;
        const Vect xc(d["cx"], d["cy"], 0.); // center
        const Scal sig = d["sig"]; // sigma
        const Scal magn = d["magn"]; // magnitude

        p.vel = [xc, sig, magn](const Vect& x) -> Vect {
          Vect dx = x - xc;
          dx[2] = 0;
          const Scal sig2 = sqr(sig);
          const Scal omz =
              magn / (2 * M_PI * sig2) * std::exp(-dx.sqrnorm() / sig2);
          return Vect(0., 0., omz);
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
  // Parses string describing a primitive.
  //   target: target name
  //   str: string to parse
  //   keys: string with keys "k0 k1 k2 ..."
  //   nreq: number of required fields
  // Format of string:
  //   name v0 v1 ...
  // where v0, v1 ... are floating point numbers,
  // If target matches `name`, returns map r={k0:v0, k1:v1, ...}
  // Otherwise, returns empty map.
  static std::map<std::string, Scal> GetMap(
      std::string target, std::string str, std::string keys, size_t nreq) {
    std::stringstream vv(str);
    vv >> std::skipws;
    std::string name;
    vv >> name;
    if (target != name) {
      return std::map<std::string, Scal>();
    }
    std::stringstream kk(keys);
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
    if (d.size() < nreq && kk) {
      throw std::runtime_error(
          "PrimList: missing field '" + k + "' in '" + str +
          "' while parsing keys '" + keys + "'");
    }
    if (d.size() < nreq && !kk) {
      throw std::runtime_error(
          "PrimList: no keys for " + std::to_string(nreq) +
          " required fields in '" + keys + "'");
    }
    if (vv && !kk) {
      throw std::runtime_error(
          "PrimList: no key for '" + std::to_string(v) + "' in '" + str +
          "' while parsing keys '" + keys + "'");
    }
    return d;
  }

 private:
  static std::string RemoveComment(std::string s) {
    return s.substr(0, s.find('#', 0));
  }

  // mods: string with allowed modifier characters (e.g. "-*")
  static std::pair<char, std::string> ExtractModifier(
      std::string s, std::string mods) {
    std::stringstream ss(s);
    ss >> std::skipws;
    std::string word;
    ss >> word;
    if (word.empty()) {
      return {' ', s};
    }
    const char c = word[0];
    if (mods.find(c) != std::string::npos) {
      std::string rest;
      std::getline(ss, rest);
      return {c, word.substr(1) + ' ' + rest};
    }
    return {' ', s};
  }

  static bool IsEmpty(std::string s) {
    s = RemoveComment(s);
    if (std::all_of(
            s.cbegin(), s.cend(), [](char c) { return std::isspace(c); })) {
      return true;
    }
    return false;
  }

  static std::string GetWord(std::string s, size_t n) {
    std::stringstream ss(s);
    ss >> std::skipws;
    std::string word;
    for (size_t i = 0; i < n + 1; ++i) {
      if (ss.good()) {
        ss >> word;
      } else {
        return "";
      }
    }
    return word;
  }

  static bool IsNumber(std::string s) {
    std::stringstream ss(s);
    Scal a;
    ss >> a;
    return !ss.fail();
  }

  static void SetDefault(
      std::map<std::string, Scal>& d, std::string key, Scal value) {
    if (!d.count(key)) {
      d[key] = value;
    }
  }
};
