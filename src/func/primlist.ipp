// Created by Petr Karnakov on 02.12.2020
// Copyright 2020 ETH Zurich

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "util/format.h"
#include "util/logger.h"

#include "primlist.h"

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
  double a;
  ss >> a;
  return !ss.fail();
}

template <class Vect_>
struct Imp {
  using Vect = Vect_;
  using Scal = typename Vect::Scal;
  static constexpr size_t dim = Vect::dim;
  using Vect2 = generic::Vect<Scal, 2>;
  using Vect3 = generic::Vect<Scal, 3>;

  using Primitive = generic::Primitive<Vect>;
  using VelocityPrimitive = generic::VelocityPrimitive<Vect>;

  static void SetDefault(
      std::map<std::string, Scal>& d, std::string key, Scal value) {
    if (!d.count(key)) {
      d[key] = value;
    }
  }

  static Scal Clamp(Scal a, Scal min, Scal max) {
    return std::max(min, std::min(max, a));
  }

  // Parses a primitive descriptor.
  // Format of the descriptor:
  //   name v0 v1 ...
  // where v0, v1 ... are floats.
  // Parameters:
  //   target_name: target name
  //   desc: primitive descriptor
  //   skeys: string with keys "k0 k1 k2 ..."
  //   nreq: number of required fields
  //   rest: buffer for values without keys
  // If target matches `name`, returns map r={k0:v0, k1:v1, ...}.
  // Fails if `keys` or `desc` contains less than `nreq` fields.
  // If target does not match `name`, returns empty map.
  // Sets `rest` to values left after parsing `keys`.
  static std::map<std::string, Scal> GetMap(
      std::string target_name, std::string desc, std::string skeys, size_t nreq,
      std::vector<Scal>* rest = nullptr) {
    std::stringstream buf_keys(skeys);
    std::stringstream buf_desc(desc);
    buf_keys >> std::skipws;
    buf_desc >> std::skipws;
    std::string name;
    buf_desc >> name;
    if (target_name != name) {
      return {};
    }

    std::vector<std::string> keys(
        std::istream_iterator<std::string>(buf_keys), {});
    std::vector<Scal> values(std::istream_iterator<Scal>(buf_desc), {});

    fassert(
        nreq <= keys.size(), //
        util::Format(
            "required {} values but only {} keys provided while parsing keys "
            "'{}'",
            nreq, keys.size(), skeys));
    fassert(
        nreq <= values.size(), //
        util::Format(
            "required {} values but only {} values provided, missing value for "
            "key '{}' in '{}' while parsing keys '{}'",
            nreq, values.size(), keys[values.size()], desc, skeys));

    std::map<std::string, Scal> res;
    size_t i;
    for (i = 0; i < std::min(keys.size(), values.size()); ++i) {
      res[keys[i]] = values[i];
    }
    if (rest) {
      (*rest) = std::vector<Scal>(values.begin() + i, values.end());
    }
    return res;
  }

  static bool ParseSphere(std::string s, size_t edim, Primitive& res) {
    res.name = "sphere";
    auto d = GetMap(res.name, s, "cx cy cz rx ry rz", 4);
    if (d.empty()) {
      d = GetMap("s", s, "cx cy cz rx ry rz", 4);
    }
    if (d.empty()) {
      return false;
    }
    SetDefault(d, "ry", d["rx"]);
    SetDefault(d, "rz", d["ry"]);

    const Vect xc(Vect3(d["cx"], d["cy"], d["cz"])); // center
    const Vect r(Vect3(d["rx"], d["ry"], d["rz"])); // radius (semi-axes)
    res.c = xc;
    res.r = r;

    res.ls = [edim, xc, r](const Vect& x) -> Scal {
      Vect xd = (x - xc) / r;
      if (edim == 2) {
        xd[2] = 0.;
      }
      return (1. - xd.sqrnorm()) * sqr(r.min());
    };

    res.inter = [xc, r](const Rect<Vect>& rect) -> bool {
      const Rect<Vect> rectbig(rect.low - r, rect.high + r);
      return rectbig.IsInside(xc);
    };
    return true;
  }

  static bool ParseRing(std::string s, Primitive& res) {
    res.name = "ring";
    auto d = GetMap(res.name, s, "cx cy cz nx ny nz r th", 8);
    if (d.empty()) {
      return false;
    }
    const Vect xc(Vect3(d["cx"], d["cy"], d["cz"])); // center
    Vect n(Vect3(d["nx"], d["ny"], d["nz"])); // normal
    n /= n.norm();
    const Scal r = d["r"]; // radius
    const Scal th = d["th"]; // thickness

    res.ls = [xc, n, r, th](const Vect& x) -> Scal {
      const Vect dx = x - xc;
      const Scal xn = dx.dot(n);
      const Scal xt = (dx - n * xn).norm();
      return sqr(th) - (sqr(xn) + sqr(xt - r));
    };
    return true;
  }

  static bool ParseBox(std::string s, size_t edim, Primitive& res) {
    res.name = "box";
    auto d = GetMap(res.name, s, "cx cy cz rx ry rz rotz", 4);
    if (d.empty()) {
      return false;
    }
    SetDefault(d, "ry", d["rx"]);
    SetDefault(d, "rz", d["ry"]);
    SetDefault(d, "rotz", 0);

    const Vect xc(Vect3(d["cx"], d["cy"], d["cz"])); // center
    const Vect r(Vect3(d["rx"], d["ry"], d["rz"])); // radius (semi-axes)
    const Scal rotz = d["rotz"]; // rotation angle around z in degrees
    const Scal rotz_cos = std::cos(-M_PI * rotz / 180.);
    const Scal rotz_sin = std::sin(-M_PI * rotz / 180.);

    res.ls = [edim, xc, r, rotz_cos, rotz_sin](const Vect& x) -> Scal {
      Vect dx = x - xc;
      const Scal x0 = dx[0] * rotz_cos - dx[1] * rotz_sin;
      const Scal x1 = dx[0] * rotz_sin + dx[1] * rotz_cos;
      dx[0] = x0;
      dx[1] = x1;
      for (size_t i = edim; i < dim; ++i) {
        dx[i] = 0;
      }
      return (r - dx.abs()).min();
    };
    return true;
  }

  static bool ParseRoundBox(std::string s, size_t edim, Primitive& res) {
    res.name = "roundbox";
    auto d = GetMap(res.name, s, "cx cy cz rx ry rz round", 7);
    if (d.empty()) {
      return false;
    }

    const Vect xc(Vect3(d["cx"], d["cy"], d["cz"])); // center
    const Vect r(Vect3(d["rx"], d["ry"], d["rz"])); // radius (semi-axes)
    const Scal round = d["round"]; // radius of rounding

    res.ls = [edim, xc, r, round](const Vect& x) -> Scal {
      Vect dx = (x - xc).abs();
      for (size_t i = edim; i < dim; ++i) {
        dx[i] = 0;
      }
      const Vect q = (dx - r + Vect(round)).max(Vect(0));
      return round - q.norm();
    };
    return true;
  }

  static bool ParseSmoothStep(std::string s, Primitive& res) {
    res.name = "smooth_step";
    // smooth step (half of diverging channel from almgren1997)
    //  +   +   +   +   +   +  +
    //              -------------       ____t
    //  +   +   +  /c   -   -  -       |
    // ------------                    |n
    //  -   -   -   -   -   -  -
    auto d = GetMap(res.name, s, "cx cy cz nx ny nz tx ty tz ln lt", 4);
    if (d.empty()) {
      return false;
    }

    const Vect xc(Vect3(d["cx"], d["cy"], d["cz"])); // center
    Vect n(Vect3(d["nx"], d["ny"], d["nz"])); // normal
    n /= n.norm();
    Vect t(Vect3(d["tx"], d["ty"], d["tz"])); // tangent
    t -= n * n.dot(t);
    t /= t.norm();
    const Scal ln = d["ln"]; // half-length in normal direction
    const Scal lt = d["lt"]; // half-length in tangential direction

    res.ls = [xc, n, t, ln, lt](const Vect& x) -> Scal {
      const Vect dx = x - xc;
      const Scal dn = n.dot(dx) / ln;
      const Scal dt = t.dot(dx) / lt;
      const Scal q = dt < -1 ? 1 : dt > 1 ? -1 : -std::sin(M_PI * 0.5 * dt);
      return (q - dn) * lt;
    };
    return true;
  }

  static bool ParseCylinder(std::string s, size_t edim, Primitive& res) {
    res.name = "cylinder";
    // c: center
    // t: axis
    // r: radius
    // t0,t1: length between center
    auto d = GetMap(res.name, s, "cx cy cz tx ty tz r t0 t1 ", 4);
    if (d.empty()) {
      return false;
    }

    if (dim == 2) {
      const Vect xc(Vect3(d["cx"], d["cy"], 0));
      const Scal r = d["r"];
      res.ls = [xc, r](const Vect& x) -> Scal { //
        return r - (x - xc).norm();
      };
    } else {
      const Vect xc(Vect3(d["cx"], d["cy"], d["cz"]));
      Vect t(Vect3(d["tx"], d["ty"], d["tz"]));
      t /= t.norm();
      const Scal r = d["r"];
      const Scal t0 = d["t0"];
      const Scal t1 = d["t1"];

      res.ls = [xc, t, r, t0, t1, edim](const Vect& x) -> Scal {
        Vect dx = x - xc;
        for (size_t i = edim; i < dim; ++i) {
          dx[i] = 0;
        }
        const Scal dt = t.dot(dx);
        const Scal dr = dx.dist(t * dt);
        Scal q = r - dr;
        if (dt < t0) q = std::min(q, dt - t0);
        if (dt > t1) q = std::min(q, t1 - dt);
        return q;
      };
    }
    return true;
  }

  template <class F>
  static Scal GetPolygonSdf(const Vect2& p, size_t npoints, F get_point) {
    // http://geomalgorithms.com/a03-_inclusion.html
    int wn = 0; // winding number
    Scal sqdist = std::numeric_limits<Scal>::max();
    size_t ibegin = 0;
    for (size_t i = 0; i + 1 < npoints; ++i) {
      const auto e = get_point(i + 1) - get_point(i);
      const auto w = p - get_point(i);
      const auto b = w - e * Clamp(w.dot(e) / e.dot(e), 0, 1);
      sqdist = std::min(sqdist, b.dot(b));
      if (get_point(i)[1] <= p[1]) {
        if (get_point(i + 1)[1] > p[1]) {
          if (e.cross_third(w) > 0) {
            ++wn;
          }
        }
      } else { // get_point(i)[1] > p[1]
        if (get_point(i + 1)[1] <= p[1]) {
          if (e.cross_third(w) < 0) {
            --wn;
          }
        }
      }
      if (get_point(i + 1) == get_point(ibegin)) { // reached end of a polygon
        ibegin = i + 2;
        ++i;
      }
    }
    return (wn == 0 ? -1 : 1) * std::sqrt(sqdist);
  }

  static void AssertPolygonPoints(const std::vector<Vect2>& points) {
    size_t ibegin = 0;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
      if (points[i + 1] == points[ibegin]) { // reached end of a polygon
        ibegin = i + 2;
        ++i;
      }
    }
    fassert_equal(
        ibegin, points.size(),
        "First and last vertex of a polygon must coincide");
  }

  static bool ParsePolygon(std::string s, size_t edim, Primitive& res) {
    res.name = "polygon";
    std::vector<Scal> coords;
    auto d = GetMap(
        res.name, s, "ox oy oz   nx ny nz   ux uy uz   n0 n1   scale", 12,
        &coords);
    if (d.empty()) {
      return false;
    }
    // coords: vertices of nonintersecting polygons stored as loops
    //         that start and end with the same point.

    fassert(
        coords.size() % 2 == 0,
        util::Format(
            "got {} coordinates, required an even number", coords.size()));
    const size_t npoints = coords.size() / 2;
    fassert(
        npoints >= 4,
        util::Format(
            "got {} two-dimensional points, at least 4 points required",
            npoints));

    std::vector<Vect2> points(npoints);
    const Scal scale = d["scale"];
    for (size_t i = 0; i < npoints; ++i) {
      points[i][0] = coords[2 * i] * scale;
      points[i][1] = coords[2 * i + 1] * scale;
    }
    AssertPolygonPoints(points);

    const Vect o(Vect3(d["ox"], d["oy"], d["oz"])); // origin
    Vect n(Vect3(d["nx"], d["ny"], d["nz"])); // normal
    n /= n.norm();
    // direction of two-dimensional x-axis
    Vect u(Vect3(d["ux"], d["uy"], d["uz"]));
    u -= n.dot(u) * n;
    u /= u.norm();
    const Vect v = n.cross(u);
    const Scal n0 = d["n0"];
    const Scal n1 = d["n1"];

    // TODO: add bounding box heuristic
    res.ls = [edim, o, n, u, v, n0, n1, points](const Vect& x) -> Scal {
      Vect dx = x - o;
      for (size_t i = edim; i < dim; ++i) {
        dx[i] = 0;
      }
      const Scal dn = n.dot(dx);
      const Scal du = u.dot(dx);
      const Scal dv = v.dot(dx);
      const Vect2 p(du, dv);
      Scal q = GetPolygonSdf(
          p, points.size(), [&points](size_t i) { return points[i]; });
      if (dn < n0) q = std::min(q, dn - n0);
      if (dn > n1) q = std::min(q, n1 - dn);
      return q;
    };
    res.inter = [](const Rect<Vect>&) -> bool { return true; };
    return true;
  }
  static bool ParsePolygon2(std::string s, Primitive& res) {
    if (dim != 2) {
      return false;
    }
    res.name = "polygon2";
    std::vector<Scal> coords;
    auto d = GetMap(res.name, s, "ox oy ux uy scale", 5, &coords);
    if (d.empty()) {
      return false;
    }
    // coords: vertices of nonintersecting polygons stored as loops
    //         that start and end with the same point.

    fassert(
        coords.size() % 2 == 0,
        util::Format(
            "got {} coordinates, required an even number", coords.size()));
    const size_t npoints = coords.size() / 2;
    fassert(
        npoints >= 4,
        util::Format(
            "got {} two-dimensional points, at least 4 points required",
            npoints));

    std::vector<Vect2> points(npoints);
    const Scal scale = d["scale"];
    for (size_t i = 0; i < npoints; ++i) {
      points[i][0] = coords[2 * i] * scale;
      points[i][1] = coords[2 * i + 1] * scale;
    }
    AssertPolygonPoints(points);

    const Vect2 o(d["ox"], d["oy"]); // origin
    // direction of local x-axis
    Vect2 u(d["ux"], d["uy"]);
    u /= u.norm();
    const Vect2 v(-u[1], u[0]);

    // TODO: add bounding box heuristic
    res.ls = [o, u, v, points](const Vect& x) -> Scal {
      Vect2 dx = Vect2(x) - o;
      const Scal du = u.dot(dx);
      const Scal dv = v.dot(dx);
      const Vect2 p(du, dv);
      Scal q = GetPolygonSdf(
          p, points.size(), [&points](size_t i) { return points[i]; });
      return q;
    };
    res.inter = [](const Rect<Vect>&) -> bool { return true; };
    return true;
  }

  static bool ParseRuled(std::string s, size_t edim, Primitive& res) {
    res.name = "ruled";
    std::vector<Scal> coords;
    auto d = GetMap(
        res.name, s, "ox oy oz   nx ny nz   ux uy uz   n0 n1   scale0 scale1",
        13, &coords);
    if (d.empty()) {
      return false;
    }
    // coords: Vertices of nonintersecting polygons stored as loops
    //         that start and end with the same point.
    //         Must be an even number of polygons,
    //         two halves used on the opposite sides of the ruled surface.

    fassert(
        coords.size() % 2 == 0,
        util::Format(
            "got {} coordinates, required an even number", coords.size()));
    const size_t npoints = coords.size() / 2;
    fassert(
        npoints >= 4,
        util::Format(
            "got {} two-dimensional points, at least 8 points required",
            npoints));
    fassert(
        npoints % 2 == 0,
        util::Format(
            "got {} two-dimensional points, required an even number", npoints));

    std::vector<Vect2> points0(npoints / 2);
    std::vector<Vect2> points1(npoints / 2);
    const Scal scale0 = d["scale0"];
    const Scal scale1 = d["scale1"];
    for (size_t i = 0; i < npoints / 2; ++i) {
      points0[i][0] = coords[2 * i] * scale0;
      points0[i][1] = coords[2 * i + 1] * scale0;
      points1[i][0] = coords[2 * (i + npoints / 2)] * scale1;
      points1[i][1] = coords[2 * (i + npoints / 2) + 1] * scale1;
    }
    for (size_t i = npoints / 2; i < npoints; ++i) {
    }
    AssertPolygonPoints(points0);
    AssertPolygonPoints(points1);

    const Vect o(Vect3(d["ox"], d["oy"], d["oz"])); // origin
    Vect n(Vect3(d["nx"], d["ny"], d["nz"])); // normal
    n /= n.norm();
    Vect u(Vect3(
        d["ux"], d["uy"], d["uz"])); // direction of two-dimensional x-axis
    u -= n.dot(u) * n;
    u /= u.norm();
    const Vect v = n.cross(u);
    const Scal n0 = d["n0"];
    const Scal n1 = d["n1"];

    // TODO: add bounding box heuristic
    res.ls = [edim, o, n, u, v, n0, n1, points0, points1](const Vect& x) {
      Vect dx = x - o;
      for (size_t i = edim; i < dim; ++i) {
        dx[i] = 0;
      }
      const Scal dn = n.dot(dx);
      const Scal du = u.dot(dx);
      const Scal dv = v.dot(dx);
      const Vect2 p(du, dv);
      const Scal alpha = Clamp((n1 - dn) / (n1 - n0), 0, 1);
      Scal q = GetPolygonSdf(
          p, points0.size(), [&points0, &points1, alpha](size_t i) {
            return alpha * points0[i] + (1 - alpha) * points1[i];
          });
      if (dn < n0) q = std::min(q, dn - n0);
      if (dn > n1) q = std::min(q, n1 - dn);
      return q;
    };
    res.inter = [](const Rect<Vect>&) -> bool { return true; };
    return true;
  }

  // Parses a list of primitives in stream buf.
  // edim: effective dimension, 2 or 3 (ignores z-component if edim=2)
  static std::vector<Primitive> GetPrimitives(std::istream& buf, size_t edim) {
    std::vector<Primitive> pp;

    buf >> std::skipws;

    // Read until eof
    while (buf) {
      std::string s;
      std::getline(buf, s);
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
            fassert(
                false,
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
      if (!r) r = ParsePolygon(s, edim, p);
      if (!r) r = ParsePolygon2(s, p);
      if (!r) r = ParseRuled(s, edim, p);

      if (p.mod_minus) {
        p.inter = [](const Rect<Vect>&) -> bool { return true; };
      }

      fassert(r, "PrimList: primitive not recognized in '" + s + "'");
      pp.push_back(p);
    }

    return pp;
  }

  static std::vector<VelocityPrimitive> GetVelocityPrimitives(
      std::istream& buf, size_t edim) {
    (void)dim;
    (void)edim;
    std::vector<VelocityPrimitive> pp;

    buf >> std::skipws;

    // Read until eof
    while (buf) {
      std::string s;
      std::getline(buf, s);
      std::map<std::string, Scal> d;

      d = GetMap("ring", s, "cx cy cz nx ny nz r th magn", 9);
      if (!d.empty()) {
        VelocityPrimitive p;
        const Vect xc(Vect3(d["cx"], d["cy"], d["cz"])); // center
        Vect n(Vect3(d["nx"], d["ny"], d["nz"])); // normal
        n /= n.norm();
        const Scal r = d["r"]; // radius
        const Scal th = d["th"]; // thickness
        const Scal magn = d["magn"]; // magnitude

        p.velocity = [xc, n, r, th, magn](const Vect& x) -> Vect {
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

      d = GetMap("gauss2d", s, "cx cy sig magn", 4);
      if (!d.empty()) {
        VelocityPrimitive p;
        const Vect xc(Vect2(d["cx"], d["cy"])); // center
        const Scal sig = d["sig"]; // sigma
        const Scal magn = d["magn"]; // magnitude

        p.velocity = [xc, sig, magn](const Vect& x) -> Vect {
          Vect dx = x - xc;
          dx[2] = 0;
          const Scal sig2 = sqr(sig);
          const Scal omz =
              magn / (2 * M_PI * sig2) * std::exp(-dx.sqrnorm() / sig2);
          if (dim == 2) {
            return Vect(Vect2(omz, 0));
          } else {
            return Vect(Vect3(0, 0, omz));
          }
        };
        pp.push_back(p);
      }
    }

    return pp;
  }
};

template <class Vect>
std::vector<generic::Primitive<Vect>> UPrimList<Vect>::GetPrimitives(
    std::istream& buf, size_t edim) {
  return Imp<Vect>::GetPrimitives(buf, edim);
}

template <class Vect>
std::vector<generic::VelocityPrimitive<Vect>>
UPrimList<Vect>::GetVelocityPrimitives(std::istream& buf, size_t edim) {
  return Imp<Vect>::GetVelocityPrimitives(buf, edim);
}
