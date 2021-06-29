// Created by Petr Karnakov on 13.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "geom/vect.h"
#include "reconst.h"

// Particle strings class.
// Single string:
// - seed particles
// - attach interface
// - compute forces from interface
// - apply constraints
// - advance
// - iteration until equilibration
// - arrays passed by pointer and length
// Multiple strings:
// - collection of independent strings
// - stored in contiguous array
// Client:
// - constructor with Par

// particle strings
// Assume 2d vectors (x[2]=0).
template <class Scal>
class PartStr {
 public:
  struct Par {
    Scal leq = 4; // length of partricle string relative to cell size
    Scal relax = 0.5; // relaxation factor
    size_t npmax = 7; // maximum number particles in string
    Scal segcirc = 0; // factor for shift to circular segment
    Scal hc = 0; // cell size [length]
    bool dn = false; // normal displacement
  };

  static constexpr size_t dim = 2;
  using Vect = generic::Vect<Scal, dim>;

  PartStr(Par par_) : par(par_) {
    Clear();
  }
  const Par& GetPar() const {
    return par;
  }
  void SetPar(Par par0) {
    par = par0;
  }
  // Number of strings
  size_t GetNumStr() const {
    return sx_.size() - 1;
  }
  // Remove all strings
  void Clear() {
    xx_.clear();
    sx_.clear();
    sx_.push_back(0);
    lx_.clear();
    ls_.clear();
    lo_.clear();
    lo_.push_back(0);
    sl_.clear();
    sl_.push_back(0);
  }
  // Adds particle string with attached interface.
  // xc: center
  // t: tangent
  // lx: nodes of all lines
  // ls: ls[l]: size of line l
  // Returns:
  // index of new string, equals GetNumStr()-1
  size_t Add(
      const Vect& xc, const Vect& t, const std::vector<Vect>& lx,
      const std::vector<size_t>& ls) {
    size_t np = par.npmax;
    // seed particles
    Scal leq = par.leq * par.hc / (np - 1); // length of segment
    for (size_t i = 0; i < np; ++i) {
      xx_.push_back(xc + t * ((Scal(i) - (np - 1) * 0.5) * leq));
    }
    sx_.push_back(xx_.size());

    size_t i = 0; // index in lx
    // loop over lines
    for (size_t l = 0; l < ls.size(); ++l) {
      // loop over nodes of line l
      for (size_t k = 0; k < ls[l]; ++k) {
        lx_.push_back(lx[i]);
        ++i;
      }
      ls_.push_back(ls[l]);
      lo_.push_back(lx_.size());
    }
    sl_.push_back(ls_.size());
    return GetNumStr() - 1;
  }
  // Computes forces ff and advances particles of single string.
  // s: string index
  // Returns:
  // r = max(ff) / (hc * relax)
  // XXX: uses static variables
  Scal Iter(size_t s) {
    assert(s < GetNumStr());
    Vect* xx = &(xx_[sx_[s]]);
    size_t sx = sx_[s + 1] - sx_[s];

    thread_local std::vector<Vect> ff;
    ff.resize(sx);

    // compute interface forces
    InterfaceForce(
        xx, sx, &(lx_[lo_[sl_[s]]]), &(ls_[sl_[s]]), sl_[s + 1] - sl_[s],
        par.segcirc, ff.data());

    ApplyConstraints(xx, sx, par.relax, par.dn, ff.data());

    // advance particles
    for (size_t i = 0; i < sx; ++i) {
      xx[i] += ff[i];
    }

    Scal r = 0.; // residual
    for (size_t i = 0; i < sx; ++i) {
      r = std::max(r, ff[i].norminf());
    }
    r /= par.hc * par.relax;
    return r;
  }
  // Iterations for single string until convergence.
  // s: string index
  // tol: tolerance
  // itermax: maximum number of iterations
  // Returns:
  // last Iter(), number of iterations
  std::pair<Scal, size_t> Run0(size_t s, Scal tol, size_t itermax) {
    Scal r = 0.;
    size_t it;
    for (it = 0; it < itermax; ++it) {
      r = Iter(s);
      if (r < tol) {
        break;
      }
    }
    return std::make_pair(r, it);
  }
  // Iterations for all strings until convergence.
  // tol: tolerance
  // itermax: maximum number of iterations
  // verb: 1: debug output for 10 strings, 2: full log
  // Returns:
  // max over all Run0()
  std::pair<Scal, size_t> Run(Scal tol, size_t itermax, int verb) {
    if (verb == 2) {
      return RunLog(tol, itermax);
    }

    Scal rm = 0.;
    size_t itm = 0;
    for (size_t s = 0; s < GetNumStr(); ++s) {
      auto rit = Run0(s, tol, itermax);
      Scal r = rit.first;
      size_t it = rit.second;
      // report
      if (verb && s % std::max<size_t>(1, GetNumStr() / 10) == 0) {
        auto fl = std::cerr.flags();
        std::cerr.precision(16);
        std::cerr << "s=" << std::setw(10) << s << " r=" << std::setw(20) << r
                  << " it=" << std::setw(5) << it << std::endl;
        std::cerr.flags(fl);
      }
      rm = std::max(rm, r);
      itm = std::max(itm, it);
    }
    return std::make_pair(rm, itm);
  }
  // Iterations for all strings until convergence, create log file.
  // tol: tolerance
  // itermax: maximum number of iterations
  // Returns:
  // max over all Run0()
  std::pair<Scal, size_t> RunLog(Scal tol, size_t itermax) {
    Scal rm = 0.;
    size_t itm = 0;
    std::ofstream f("part.log");
    f.precision(25);
    f << "string iter k res" << std::endl;
    for (size_t s = 0; s < GetNumStr(); ++s) {
      Scal r = 0;
      size_t it = 0;
      for (it = 0; it < itermax; ++it) {
        r = Iter(s);

        // log
        f << s << ' ' << it << ' ' << GetCurv(s) << ' ' << r << '\n';

        if (r < tol) {
          break;
        }
      }

      rm = std::max(rm, r);
      itm = std::max(itm, it);
    }
    return std::make_pair(rm, itm);
  }
  // Curvature of single string
  // s: string index
  Scal GetCurv(size_t s) const {
    const Vect* xx = &(xx_[sx_[s]]);
    size_t sx = sx_[s + 1] - sx_[s];
    return PartK(xx, sx);
  }
  // Returns single string.
  // s: string index
  // Output:
  // array of positions, length
  std::pair<const Vect*, size_t> GetStr(size_t s) {
    return std::make_pair(&(xx_[sx_[s]]), sx_[s + 1] - sx_[s]);
  }
  struct Inter {
    const Vect* x; // nodes
    const size_t* s; // sizes
    size_t n; // number of lines
  };
  // Returns interface for single string.
  // s: string index
  // Output:
  // Inter
  Inter GetInter(size_t s) {
    return {&(lx_[lo_[sl_[s]]]), &(ls_[sl_[s]]), sl_[s + 1] - sl_[s]};
  }

 private:
  Par par;

  // size of arrays sx_, sl_, sk_ is the number of strings
  // positions of string i start from xx_[sx_[i]]
  // lines attached to string i start from ll_[sl_[i]]

  std::vector<Vect> xx_; // particle positions
  std::vector<size_t> sx_; // xx index, ending with sx.size()
  std::vector<Vect> lx_; // nodes of all lines
  std::vector<size_t> ls_; // ls_[l]: size of line l
  std::vector<size_t> lo_; // lo_[l]: first node of line l, index in lx_
                           // lo_[-1]: lx_.size()
  std::vector<size_t> sl_; // sl_[s]: first line of string s, index in ls_,lo_
                           // sl_[-1]: ls_.size()

  using R = Reconst<Scal>;

 public:
  // Curvature of particle string.
  // xx: array of positions
  // sx: size of xx
  static Scal PartK(const Vect* xx, size_t sx) {
    assert(sx % 2 == 1 && sx >= 3);
    size_t i = (sx - 1) / 2;
    Vect x = xx[i];
    Vect xm = xx[i - 1];
    Vect xp = xx[i + 1];
    Vect dm = x - xm;
    Vect dp = xp - x;
    Scal lm = dm.norm();
    Scal lp = dp.norm();
    Scal lmp = lm * lp;
    Scal k = std::sqrt(2. * (lmp - dm.dot(dp))) / lmp;
    if (dm.cross_third(dp) > 0.) {
      k = -k;
    }
    return k;
  }
  // Circular profile along segment.
  // k: curvature
  // l: segment half-length
  // d: distance from segment center to target point (d<l)
  static Scal SegCirc(Scal k, Scal l, Scal d) {
    d = std::min(d, l); // XXX: adhoc
    k = std::min(k, 1. / l);
    Scal a1 = 1. - sqr(k * l);
    Scal a2 = 1. - sqr(k * d);
    Scal t1 = std::sqrt(a1);
    Scal t2 = std::sqrt(a2);
    return k * (sqr(l) - sqr(d)) / (t1 + t2);
  };
  // Force on particles from interface.
  // xx: array of positions
  // sx: size of xx
  // lx: nodes of lines
  // ls: sizes of lines
  // sl: size of ls
  // segcirc: factor for shift to circular segment
  // anglim: limit for angle between string and interface
  // ft: force type
  // Output:
  // ff: forces
  static void InterfaceForce(
      const Vect* xx, size_t sx, const Vect* lx, const size_t* ls, size_t sl,
      Scal segcirc, Vect* ff) {
    if (sx == 0) {
      return;
    }
    assert(sx % 2 == 1);

    // clear force
    for (size_t i = 0; i < sx; ++i) {
      ff[i] = Vect(0);
    }

    if (sl == 0) {
      return;
    }

    // current curvature
    Scal crv = PartK(xx, sx);

    // find nearest segment over all interface lines
    // x: target point
    // Returns:
    // e: segment ends
    auto findnear = [&lx, &ls, &sl](const Vect& x) -> std::array<Vect, 2> {
      std::array<Vect, 2> en = {Vect(0), Vect(0)}; // result
      Scal dn = std::numeric_limits<Scal>::max(); // sqrdist to nearest

      // loop over interface lines
      size_t b = 0; // index in lx
      for (size_t l = 0; l < sl; ++l) {
        // loop over segments of line l
        for (size_t k = 0; k < ls[l]; ++k) {
          size_t kp = (k + 1 == ls[l] ? 0 : k + 1);
          Scal d = x.sqrdist(R::GetNearest(x, lx[b + k], lx[b + kp]));
          if (d < dn) {
            dn = d;
            en = {lx[b + k], lx[b + kp]};
          }
          if (ls[l] == 2) {
            break;
            ;
          }
        }
        b += ls[l];
      }
      return en;
    };

    // Apply shift from segment to circle
    // e: segment ends
    // x: curent point on segment
    auto shsegcirc = [crv, segcirc](const std::array<Vect, 2>& e, Vect& x) {
      if (segcirc != 0.) {
        // outer normal
        Vect n = e[1] - e[0];
        n = Vect(n[1], -n[0]);
        Scal nn = n.norm();
        if (nn < 1e-10) { // XXX: adhoc
          return;
        }
        n /= nn;
        // center
        Vect xc = (e[0] + e[1]) * 0.5;
        // distance from center
        Scal dc = xc.dist(x);
        // max distance from center
        Scal mdc = e[0].dist(xc);
        // shift from line to circle
        Scal s = SegCirc(crv, mdc, dc);
        s *= segcirc;
        x += n * s;
      }
    };

    for (size_t i = 0; i < sx; ++i) {
      const Vect& x = xx[i];
      auto e = findnear(x);
      Vect xn = R::GetNearest(x, e[0], e[1]);
      shsegcirc(e, xn);
      ff[i] += xn - x;
    }
  }
  // Oriented angle from (x1-x0) to (x2-x1)
  static Scal GetAn(Vect x0, Vect x1, Vect x2) {
    Vect dm = x1 - x0;
    Vect dp = x2 - x1;
    Scal lm = dm.norm();
    Scal lp = dp.norm();
    Scal sin = dm.cross_third(dp) / (lm * lp);
    return std::asin(sin);
  };
  // Oriented angle from (x[i]-x[i-1]) to (x[i+1]-x[i])
  static Scal GetAn(const Vect* x, size_t i) {
    return GetAn(x[i - 1], x[i], x[i + 1]);
  };
  // Apply constraints on force.
  // xx: array of positions
  // sx: size of xx
  // relax: relaxation factor for force
  // dn: 1: enable normal displacement, 0: fix central particle
  static void ApplyConstraints(
      const Vect* xx, size_t sx, Scal relax, bool dn, Vect* ff) {
    // Rotates vector by pi/2
    // x: vector of plane coordinates
    auto rr = [](const Vect& x) { return Vect(-x[1], x[0]); };
    // Rotates vector to angle 'a' given by unit vector
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto re = [](const Vect& x, const Vect& e) {
      return Vect(x[0] * e[0] - x[1] * e[1], x[0] * e[1] + x[1] * e[0]);
    };
    // Rotates vector to angle '-a' with 'a' given by unit vector
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto rem = [](const Vect& x, const Vect& e) {
      return Vect(x[0] * e[0] + x[1] * e[1], -x[0] * e[1] + x[1] * e[0]);
    };
    // Returns vector at angle a
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto ra = [](Scal a) { return Vect(std::cos(a), std::sin(a)); };

    // relaxation of force
    for (size_t i = 0; i < sx; ++i) {
      ff[i] *= relax;
    }

    // new positions
    thread_local std::vector<Vect> xxn; // dx/dtheta
    xxn.resize(sx);
    // initialize
    for (size_t i = 0; i < sx; ++i) {
      xxn[i] = xx[i];
    }

    // Note:
    // New positions without constraints would be xxn=xx+ff.
    // ff acts as correction of positions.
    // Constraints are applied separately:
    // 1) displacement of center:
    //    xxn += dx
    //    ff -= dx
    // 2) alpha
    //    xxn = X(alpha)
    //    ff -= dx
    // 3) theta
    //    xxn = X(theta)
    //    ff -= dx

    // displacement of center
    if (dn) { // average over all particles
      Vect dx(0);
      for (size_t i = 0; i < sx; ++i) {
        dx += ff[i];
      }
      dx /= sx;
      dx[0] = 0.;

      // apply
      for (size_t i = 0; i < sx; ++i) {
        xxn[i] += dx;
        ff[i] -= dx;
      }
    } else { // by central particle
      Vect dx(0);
      dx[1] += ff[(sx - 1) / 2][1];

      // apply
      for (size_t i = 0; i < sx; ++i) {
        xxn[i] += dx;
        ff[i] -= dx;
      }
    }

    // alpha: angle between x-axis and normal
    {
      // derivative of positions by alpha
      thread_local std::vector<Vect> xa; // dx/dalpha
      xa.resize(sx);
      // central
      const size_t ic = (sx - 1) / 2;
      xa[ic] = Vect(0);
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        Vect d;
        // forward
        i = ic + q;
        d = xxn[i] - xxn[i - 1];
        xa[i] = xa[i - 1] + rr(d);
        // backward
        i = ic - q;
        d = xxn[i + 1] - xxn[i];
        xa[i] = xa[i + 1] - rr(d);
      }

      // projection of f
      Scal fa = 0.; // f dot xa
      Scal na = 0.; // xa dot xa
      for (size_t i = 0; i < sx; ++i) {
        fa += ff[i].dot(xa[i]);
        na += xa[i].dot(xa[i]);
      }

      // correction of alpha
      Scal da = fa / na;

      // vector at angle da
      Vect ea = ra(da);

      // segment vectors
      thread_local std::vector<Vect> dd;
      dd.resize(sx);
      dd[ic] = Vect(0);
      // initialize from xxn
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = xxn[i] - xxn[i - 1];
        // backward
        i = ic - q;
        dd[i] = xxn[i + 1] - xxn[i];
      }

      // apply da
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = re(dd[i], ea);
        // backward
        i = ic - q;
        dd[i] = re(dd[i], ea);
      }

      // restore new xxn from segments
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        Vect dx;
        // forward
        i = ic + q;
        dx = xxn[i - 1] + dd[i] - xxn[i];
        xxn[i] += dx;
        ff[i] -= dx;
        // backward
        i = ic - q;
        dx = xxn[i + 1] - dd[i] - xxn[i];
        xxn[i] += dx;
        ff[i] -= dx;
      }
    }

    // theta: angle between segments
    {
      // derivative of positions by theta
      thread_local std::vector<Vect> xt; // dx/dtheta
      xt.resize(sx);
      // central
      const size_t ic = (sx - 1) / 2;
      xt[ic] = Vect(0);
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        Vect d;
        // forward
        i = ic + q;
        d = xxn[i] - xxn[i - 1];
        xt[i] = xt[i - 1] + rr(d) * (q - 0.5);
        // backward
        i = ic - q;
        d = xxn[i + 1] - xxn[i];
        xt[i] = xt[i + 1] + rr(d) * (q - 0.5);
      }

      // projection of f
      Scal ft = 0.; // f dot xt
      Scal nt = 0.; // xt dot xt
      for (size_t i = 0; i < sx; ++i) {
        ft += ff[i].dot(xt[i]);
        nt += xt[i].dot(xt[i]);
      }

      // correction of theta
      Scal dt = ft / nt;
      // vector at angle dt/2
      Vect eth = ra(dt);
      // vector at angle dt
      Vect et = re(eth, eth);

      // segment vectors
      thread_local std::vector<Vect> dd;
      dd.resize(sx);
      dd[ic] = Vect(0);
      // initialize from xxn
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = xxn[i] - xxn[i - 1];
        // backward
        i = ic - q;
        dd[i] = xxn[i + 1] - xxn[i];
      }

      // apply dt
      {
        Vect e = eth;
        for (size_t q = 1; q <= ic; ++q) {
          size_t i;
          // forward
          i = ic + q;
          dd[i] = re(dd[i], e);
          // backward
          i = ic - q;
          dd[i] = rem(dd[i], e);
          // next
          e = re(e, et);
        }
      }

      // restore new xxn from segments
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        Vect dx;
        // forward
        i = ic + q;
        dx = xxn[i - 1] + dd[i] - xxn[i];
        xxn[i] += dx;
        ff[i] -= dx;
        // backward
        i = ic - q;
        dx = xxn[i + 1] - dd[i] - xxn[i];
        xxn[i] += dx;
        ff[i] -= dx;
      }
    }

    // convert to position correction
    for (size_t i = 0; i < sx; ++i) {
      ff[i] = xxn[i] - xx[i];
    }
  }
};
