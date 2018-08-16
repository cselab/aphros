#pragma once

#include <exception>
#include <array>
#include <memory>
#include <limits>
#include <vector>
#include <iomanip>

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
    Scal leq = 4.; // equilibrium length of string relative to cell size
    Scal relax = 0.5; // relaxation factor
    int constr = 1;  // constraint type
    size_t npmax = 11; // maximum number particles in string
    Scal segcirc = 1.; // factor for shift to circular segment
    Scal hc; // cell size [length]
    // specific for Constr0
    Scal kstr = 1.;  // factor for stretching 
    Scal kbend = 1.;  // factor for bending 
    Scal kattr = 1.;  // factor for attraction 
    bool bendmean = true; // 1: bending to mean angle, 0: to straight angle
    // specific for Constr1
    Scal ka = 1.;  // factor for angle between normal and x-axis 
    Scal kt = 1.;  // factor for angle between segments
    Scal kx = 1.;  // factor for position
    Scal tmax = 180.; // limit for total theta of string [deg]
    Scal dtmax = 10.; // limit for dt [deg]
    Scal anglim = 90.; // limit for dt [deg]
  };

  static constexpr size_t dim = 3;
  using Vect = GVect<Scal, dim>;

  PartStr(std::shared_ptr<Par> par) : par(par) {
    Clear();
  }
  Par* GetPar() { return par.get(); }
  // Number of strings
  size_t GetNumStr() const {
    return sx_.size() - 1;
  }
  // Remove all strings
  void Clear() {
    xx_.clear();
    sx_.clear();
    sx_.push_back(0);
    ll_.clear();
    sl_.clear();
    sl_.push_back(0);
  }
  // Adds particle string with attached interface.
  // xc: center 
  // t: tangent
  // ll: array of lines
  // sl: size of ll
  // Returns:
  // index of new string, equals GetNumStr()-1
  size_t Add(const Vect& xc, const Vect& t, 
             const std::array<Vect, 2>* ll, size_t sl) {
    size_t np = par->npmax;
    Scal leq = par->leq * par->hc / np; // length of segment
    for (size_t i = 0; i < np; ++i) {
      xx_.push_back(xc + t * ((Scal(i) - (np - 1) * 0.5) * leq));
    }
    for (size_t i = 0; i < sl; ++i) {
      ll_.push_back(ll[i]);
    }
    sx_.push_back(xx_.size());
    sl_.push_back(ll_.size());
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

    static std::vector<Vect> ff;
    ff.resize(sx);
    
    // compute interface forces
    InterfaceForce(xx, sx, &(ll_[sl_[s]]), sl_[s + 1] - sl_[s],
                   par->segcirc, par->anglim * M_PI / 180., ff.data());

    // apply constraints
    int cs = par->constr;
    if (cs == 0) {
      Constr0(xx, sx, par->kstr, par->leq * par->hc / par->npmax,
              par->kbend, par->bendmean, par->relax, ff.data());
    } else if (cs == 1) {
      Constr1(xx, sx, par->kattr, 
              par->kbend, par->kstr, par->relax, 
              par->tmax * M_PI / 180., par->dtmax * M_PI / 180.,
              ff.data());
    } else {
      throw std::runtime_error("Unknown constr=" + std::to_string(cs));
    }

    // advance particles
    for (size_t i = 0; i < sx; ++i) {
      xx[i] += ff[i];
    }

    Scal r = 0.; // residual
    for (size_t i = 0; i < sx; ++i) {
      r = std::max(r, ff[i].norminf());
    }
    r /= par->hc * par->relax;
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
  // verb: 1: debug output
  // Returns:
  // max over all Run0()
  std::pair<Scal, size_t> Run(Scal tol, size_t itermax, bool verb) {
    Scal rm = 0.;
    size_t itm = 0;
    for (size_t s = 0; s < GetNumStr(); ++s) {
      auto rit = Run0(s, tol, itermax);
      Scal r = rit.first;
      size_t it = rit.second;
      // report
      if (verb && s % std::max<size_t>(1, GetNumStr() / 10) == 0) {
        auto fl = std::cout.flags();
        std::cout.precision(16);
        std::cout 
            << "s=" << std::setw(10) << s 
            << " r=" << std::setw(20) << r 
            << " it=" << std::setw(5) << it 
            << std::endl;
        std::cout.flags(fl);
      }
      rm = std::max(rm, r);
      itm = std::max(itm, it);
    }
    return std::make_pair(rm, itm);
  }
  // Curvature of single string
  // s: string index
  Scal GetCurv(size_t s) {
    Vect* xx = &(xx_[sx_[s]]);
    size_t sx = sx_[s + 1] - sx_[s];
    return PartK(xx, sx);
  }
  // Returns single string.
  // s: string index
  // Output:
  // array of positions, length
  std::pair<Vect*, size_t> GetStr(size_t s) {
    return std::make_pair(&(xx_[sx_[s]]), sx_[s + 1] - sx_[s]);
  }

 private:
  std::shared_ptr<Par> par;

  // size of arrays sx_, sl_, sk_ is the number of strings
  // positions of string i start from xx_[sx_[i]]
  // lines attached to string i start from ll_[sl_[i]]

  std::vector<Vect> xx_; // particle positions
  std::vector<size_t> sx_; // xx index, ending with sx.size() 
  std::vector<std::array<Vect, 2>> ll_; // interface lines
  std::vector<size_t> sl_; // sl index, ending with sl.size() 

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
    Scal t1 = std::sqrt(1. - sqr(k) * sqr(l));
    Scal t2 = std::sqrt(1. - sqr(k) * sqr(d));
    return k * (sqr(l) - sqr(d)) / (t1 + t2);
  };
  // Force on particles from interface.
  // xx: array of positions
  // sx: size of xx
  // ll: array of lines
  // sl: size of ll
  // segcirc: factor for shift to circular segment 
  // anglim: limit for angle between string and interface
  // Output:
  // ff: forces
  static void InterfaceForce(const Vect* xx, size_t sx, 
                      const std::array<Vect, 2>* ll, size_t sl, 
                      Scal segcirc, Scal anglim, Vect* ff) {
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

    // find nearest segment
    // x: target point
    // Returns:
    // jn: index of nearest segment ll[jn]
    auto findnear = [sl,ll](const Vect& x) -> size_t {
      Scal dn = std::numeric_limits<Scal>::max();  // sqrdist to nearest
      size_t jn = 0;
      // find distance to nearest segment
      for (size_t j = 0; j < sl; ++j) {
        auto& e = ll[j];
        Scal dl = x.sqrdist(R::GetNearest(x, e[0], e[1]));
        if (dl < dn) {
          dn = dl;
          jn = j;
        }
      }
      return jn;
    };

    /*
    // find nearest segment with angle limit
    // x: target point
    // n: normal to string
    // ang: angle limit (skip if exceeded)
    // Returns:
    // jn: index of nearest segment ll[jn]
    auto findnearang = [sl,ll](const Vect& x, const Vect& n, Scal ang) {
      Scal dn = std::numeric_limits<Scal>::max();  // sqrdist to nearest
      size_t jn = size_t(-1);
      // find distance to nearest segment
      for (size_t j = 0; j < sl; ++j) {
        auto& e = ll[j];
        Vect ne = e[1] - e[0];
        ne = Vect(ne[1], -ne[0], 0.);
        if (ne.dot(n) / (ne.norm() * n.norm()) >= std::cos(ang)) {
          Scal dl = x.sqrdist(R::GetNearest(x, e[0], e[1]));
          if (dl < dn) {
            dn = dl;
            jn = j;
          }
        }
      }
      return jn;
    };
    */

    /*
    // find segment with nearest angle
    // x: point
    // t: direction
    // Returns:
    // jn: index of segment ll[jn] with one end e minimizing angle
    //     between (x-e) and t (maximizing t.dot(x-e) / (x-a).norm())
    auto findang = [sl,ll](const Vect& x, const Vect& t, Scal sgn) -> size_t {
      Scal sn = -std::numeric_limits<Scal>::max();  // scalar product with t
      size_t jn = -1;
      for (size_t j = 0; j < sl; ++j) {
        for (size_t k = 0; k < 2; ++k) {
          auto e = ll[j][k];
          auto s = t.dot(e - x) / (e - x).norm();
          if (s > sn && t.cross_third(e - x) * sgn >= 0.) {
            sn = s;
            jn = j;
          }
        }
      }
      assert(jn != -1);
      return jn;
    };
    */

    /*
    // find nearest segment with distance along line.
    // x: point
    // t: direction
    // Returns:
    // jn: index of segment ll[jn] such that line x+t*a intersects the segment
    //     and the intersection point provides minimal distance to x
    // jn=-1 if no intersecting segments
    auto findint = [sl,ll,leq](const Vect& x, const Vect& t) -> size_t {
      Scal dn = std::numeric_limits<Scal>::max();  // sqrdist to nearest
      size_t jn = -1;
      // find distance to nearest segment
      for (size_t j = 0; j < sl; ++j) {
        auto& e = ll[j];
        Vect xi;
        bool in = R::GetInterLine(e, x, t, xi);
        Scal d = x.sqrdist(xi);
        if (in && d < dn && d < sqr(4 * leq)) {
          dn = d;
          jn = j;
        }
      }
      return jn;
    };
    */

    // Apply shift from segment to circle
    // j: segment idx
    // x: curent point on segment
    auto shsegcirc = [ll,sl,crv,segcirc](size_t j, Vect& x) {
      if (segcirc != 0.) {
        auto& e = ll[j];
        // outer normal
        Vect n = e[1] - e[0];
        n = Vect(n[1], -n[0], 0.);
        n /= n.norm();
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

    /*
    // central: 
    // attraction to nearest interface
    const size_t ic = (sx - 1) / 2;
    {
      const Vect& x = xx[ic];
      size_t jn = findnear(x);
      auto& e = ll[jn];
      Vect xn = R::GetNearest(x, e[0], e[1]);
      shsegcirc(jn, xn);
      ff[ic] = xn - x;
    }

    // forward and backward:
    for (size_t q = 1; q <= ic; ++q) {
      for (int g : {-1, 1}) {
        size_t i = ic + q * g;
        Vect x = xx[i]; // target point
        Vect n = xx[i] - xx[i - g]; // normal to string
        n = Vect(-n[1], n[0], 0.) * g;
        size_t jn = findnearang(xx[i], n, anglim); // nearest segment
        if (jn != size_t(-1)) {
          auto& e = ll[jn];
          Vect xn = R::GetNearest(x, e[0], e[1]);
          shsegcirc(jn, xn);
          ff[i] += xn - x;  
        }
      }
    }
    */

    // particle attracted to nearest interface
    for (size_t i = 0; i < sx; ++i) {
      const Vect& x = xx[i];
      size_t jn = findnear(x);
      auto& e = ll[jn];
      Vect xn = R::GetNearest(x, e[0], e[1]);
      shsegcirc(jn, xn);
      ff[i] = (xn - x);
    }

    #if 0
    // segment turned to match angle of nearest interface
    for (size_t im = 0; im + 1 < sx; ++im) {
      size_t ip = im + 1;
      const Vect& xm = xx[im];
      const Vect& xp = xx[ip];
      Vect x = (xm + xp) * 0.5; // segment center
      size_t jn = findnear(x);
      auto& e = ll[jn];
      Vect ne(e[0][1] - e[1][1], e[1][0] - e[0][0], 0.); // normal to interface
      Vect ns(xp[1] - xm[1], xm[0] - xp[0], 0.); // normal to segment
      if (ne.dot(ns) > 0.) {
        ne /= ne.norm();
        ns /= ns.norm();
        Scal sin = ne.cross_third(ns); // sin of angle between normals
        (void)sin;
        //Vect xn = R::GetNearest(x, e[0], e[1]);
        // torque
        Scal t = segcirc * sin;
        // segment length
        Scal l = xp.dist(xm);
        // forces
        Vect fm = ns * (t / (2. * l));
        Vect fp = -ns * (t / (2. * l));
        ff[im] += fm;
        ff[ip] += fp;
      }
    }
    #endif

    (void)anglim;
    (void)crv;
    (void)segcirc;
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
    return GetAn(x[i-1], x[i], x[i+1]);
  };
  // Stretching and bending forces.
  // xx: array of positions
  // sx: size of xx
  // kstr: factor for stretching
  // leq: equilibrium length
  // kbend: factor for bending
  // bendmean: 1: bending to mean angle, 0: bending to straight angle
  // relax: relaxation factor
  // Output:
  // ff: appended with forces
  static void Constr0(const Vect* xx, size_t sx, 
                      Scal kstr, Scal leq, Scal kbend, 
                      bool bendmean, Scal relax, Vect* ff) {
    // stretching 
    for (size_t i = 0; i < sx - 1; ++i) {
      Vect dx = xx[i + 1] - xx[i];
      Vect f = dx * (kstr * (1. - leq / dx.norm()));
      ff[i] += f;
      ff[i + 1] -= f;
    }

    // mean angle
    Scal anm = (GetAn(xx,1) + GetAn(xx,2) + GetAn(xx,3)) / 3.; 

    // bending 
    for (size_t i = 1; i < sx - 1; ++i) {
      size_t im = i - 1;
      size_t ip = i + 1;
      Vect x = xx[i];
      Vect xm = xx[im];
      Vect xp = xx[ip];
      Vect dm = x - xm;
      Vect dp = xp - x;
      Scal lm = dm.norm();
      Scal lp = dp.norm();

      // normals to segments, <nm,dm> positively oriented
      Vect nm = Vect(dm[1], -dm[0], 0.) / lm; 
      Vect np = Vect(dp[1], -dp[0], 0.) / lp;
      // torque
      Scal t = kbend * lm * lp * (GetAn(xx,i) - anm * (bendmean ? 1. : 0.));
      // forces
      Vect fm = nm * (t / (2. * lm));
      Vect fp = np * (t / (2. * lp));
      // apply
      ff[im] += fm;
      ff[ip] += fp;
      ff[i] -= (fm + fp);
    }

    // relaxation
    for (size_t i = 0; i < sx; ++i) {
      ff[i] *= relax;
    }

    // freeze central particle
    Vect fc = ff[(sx - 1) / 2];
    for (size_t i = 0; i < sx; ++i) {
      ff[i] -= fc;
    }
  }
  // Apply exact constraints on force
  // xx: array of positions
  // sx: size of xx
  // ka: relaxation factor for angle between normal and x-axis
  // kt: relaxation factor for angle between segments
  // kx: relaxation factor for position 
  // relax: relaxation factor for force
  // XXX: uses static variables
  static void Constr1(const Vect* xx, size_t sx, 
                      Scal ka, Scal kt, Scal kx, Scal relax, 
                      Scal tmax, Scal dtmax,
                      Vect* ff) {
    // Rotates vector by pi/2
    // x: vector of plane coordinates
    auto rr = [](const Vect& x) {
      return Vect(-x[1], x[0], 0.);
    };
    // Rotates vector to angle 'a' given by unit vector
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto re = [](const Vect& x, const Vect& e) {
      return Vect(x[0] * e[0] - x[1] * e[1], x[0] * e[1] + x[1] * e[0], 0.);
    };
    // Rotates vector to angle '-a' with 'a' given by unit vector
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto rem = [](const Vect& x, const Vect& e) {
      return Vect(x[0] * e[0] + x[1] * e[1], -x[0] * e[1] + x[1] * e[0], 0.);
    };
    // Returns vector at angle a
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto ra = [](Scal a) {
      return Vect(std::cos(a), std::sin(a), 0.);
    };

    // relaxation of force
    for (size_t i = 0; i < sx; ++i) {
      ff[i] *= relax;
    }

    // new positions
    static std::vector<Vect> xxn; // dx/dtheta
    xxn.resize(sx);
    // initialize
    for (size_t i = 0; i < sx; ++i) {
      xxn[i] = xx[i];
    }

    // Note:
    // New positions without contraints would be xxn=xx+ff.
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
    {
      Vect dx(0);
      for (size_t i = 0; i < sx; ++i) {
        dx += ff[i];
      }
      dx *= kx / sx;
      dx[0] = 0.;

      // apply
      for (size_t i = 0; i < sx; ++i) {
        xxn[i] += dx;
        ff[i] -= dx;
      }
    }

    // alpha: angle between x-axis and normal
    {
      // derivative of positions by alpha
      static std::vector<Vect> xa; // dx/dalpha
      xa.resize(sx);
      // central 
      const size_t ic = (sx - 1) / 2; 
      xa[ic] = Vect(0.);
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
      Scal fa = 0.;  // f dot xa
      Scal na = 0.; // xa dot xa
      for (size_t i = 0; i < sx; ++i) {
        fa += ff[i].dot(xa[i]); 
        na += xa[i].dot(xa[i]);
      }

      // correction of alpha
      Scal da = fa / na * ka;

      // vector at angle da
      Vect ea = ra(da);

      // segment vectors 
      static std::vector<Vect> dd;
      dd.resize(sx);
      dd[ic] = Vect(0.);
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
      static std::vector<Vect> xt; // dx/dtheta
      xt.resize(sx);
      // central 
      const size_t ic = (sx - 1) / 2; 
      xt[ic] = Vect(0.);
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
      Scal ft = 0.;  // f dot xt
      Scal nt = 0.; // xt dot xt
      for (size_t i = 0; i < sx; ++i) {
        ft += ff[i].dot(xt[i]); 
        nt += xt[i].dot(xt[i]);
      }

      // correction of theta
      Scal dt = ft / nt * kt;
      // limabs
      // la > 0.
      /*
      auto limabs = [](Scal a, Scal la) {
        return std::min(std::abs(a), la) * (a > 0. ? 1. : -1);
      };
      dt = limabs(dt, dtmax / (sx - 1));
      Scal t = GetAn(xx, (sx - 1) / 2);
      dt = limabs(t + dt, tmax / (sx - 1)) - t;
      */
      (void) tmax;
      (void) dtmax;

      // vector at angle dt/2
      Vect eth = ra(dt);
      // vector at angle dt
      Vect et = re(eth, eth);

      // segment vectors 
      static std::vector<Vect> dd;
      dd.resize(sx);
      dd[ic] = Vect(0.);
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
