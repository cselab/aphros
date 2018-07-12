#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>
#include <vector>

#include "geom/block.h"
#include "dump/dumper.h"
#include "geom/vect.h"
#include "reconst.h"

// attraction to exact sphere
#define ADHOC_ATTR 0
// normal from exact sphere
#define ADHOC_NORM 0

// particle strings
template <class Scal>
class GPartStr {
 public:
  static constexpr size_t dim = 3;
  using Vect = GVect<Scal, dim>;
  using R=GReconst<Scal>;

  // Curvature of particle string.
  // xx: array of positions
  // sx: size of xx
  static Scal PartK(const Vect* xx, size_t sx) {
    assert(sx % 2 == 1);
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
  // Output:
  // ff: forces
  static void InterfaceForce(const Vect* xx, size_t sx, 
                      const std::array<Vect, 2>* ll, size_t sl, 
                      Scal segcirc, Vect* ff) {
    if (!sx) {
      return;
    }
    assert(sx % 2 == 1);

    // clear force
    for (size_t i = 0; i < sx; ++i) {
      ff[i] = Vect(0);
    }

    // current curvature
    Scal k = PartK(xx, sx);

    // attracting spring to nearest point on nearest line
    for (size_t i = 0; i < sx; ++i) {
      Vect x = xx[i];

      if (sl) {
        Vect xn;  // nearest point
        Scal dn = std::numeric_limits<Scal>::max();  // sqrdist to nearest
        size_t jn = 0;
        for (size_t j = 0; j < sl; ++j) {
          auto& e = ll[j];
          Vect xl = R::GetNearest(x, e[0], e[1]);
          Scal dl = xl.sqrdist(x);

          if (dl < dn) {
            xn = xl;
            dn = dl;
            jn = j;
          }
        }

        #if ADHOC_ATTR
        auto bc = ll[sl-1][0];
        auto br = ll[sl-1][1][0];
        auto dx = x - bc;
        dx /= dx.norm();
        xn = bc + dx * br;
        #endif 

        if (segcirc != 0.) {
          auto& e = ll[jn];
          // outer normal
          Vect n = e[1] - e[0];
          n = Vect(n[1], -n[0], 0.);
          n /= n.norm();
          // center
          Vect xc = (e[0] + e[1]) * 0.5;
          // distance from center
          Scal dc = xc.dist(xn);
          // max distance from center
          Scal mdc = e[0].dist(xc);

          // shift from line to circle
          Scal s = SegCirc(k, mdc, dc);  
          s *= segcirc;
          xn += n * s;
        }

        ff[i] += xn - x; 
      }
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
  // hm: cell size
  // relax: relaxation factor for force
  // XXX: uses static variables
  static void Constr1(const Vect* xx, size_t sx, 
                      Scal ka, Scal kt, Scal kx, Scal hm, Scal relax,
                      Vect* ff) {
    // relaxation
    for (size_t i = 0; i < sx; ++i) {
      ff[i] *= relax;
    }

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

    // alpha: angle between x-axis and normal
    // theta: angle between segments
    
    // derivatives of positions by angles
    static std::vector<Vect> xa(sx); // dx/dalpha
    static std::vector<Vect> xt(sx); // dx/dtheta
    // central 
    const size_t ic = (sx - 1) / 2; 
    xa[ic] = Vect(0.);
    xt[ic] = Vect(0.);
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      Vect d;
      // forward 
      i = ic + q;
      d = xx[i] - xx[i - 1];
      xa[i] = xa[i - 1] + rr(d);
      xt[i] = xt[i - 1] + rr(d) * (q - 0.5);
      // backward
      i = ic - q;
      d = xx[i + 1] - xx[i];
      xa[i] = xa[i + 1] - rr(d);
      xt[i] = xt[i + 1] + rr(d) * (q - 0.5);
    }

    // correction of angles
    Scal da = 0.;
    Scal dt = 0.;
    for (size_t i = 0; i < sx; ++i) {
      da += ff[i].dot(xa[i]); // scale hm*hm
      dt += ff[i].dot(xt[i]);
    }

    // rescale to 1
    da /= hm * hm; 
    dt /= hm * hm; 
    
    // relaxation
    da *= ka;
    dt *= kt;

    // vector at angle da
    Vect ea = ra(da);
    // vector at angle dt/2
    Vect eth = ra(dt);
    // vector at angle dt
    Vect et = re(eth, eth);

    // segment vectors 
    static std::vector<Vect> dd(sx);
    dd[ic] = Vect(0.);
    // initialize from xx
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      // forward
      i = ic + q;
      dd[i] = xx[i] - xx[i - 1];
      // backward
      i = ic - q;
      dd[i] = xx[i + 1] - xx[i];
    }

    // apply da
    {
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = re(dd[i], ea);
        // backward
        i = ic - q;
        dd[i] = re(dd[i], ea);
      }
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

    // displacement of center
    Vect dx(0);
    for (size_t i = 0; i < sx; ++i) {
      dx += ff[i];
    }
    dx *= kx;

    // restore new xx from segments, store in ff
    ff[ic] = xx[ic];
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      // forward
      i = ic + q;
      ff[i] = ff[i - 1] + dd[i];
      // backward
      i = ic - q;
      ff[i] = ff[i + 1] - dd[i];
    }

    // apply dx
    for (size_t i = 0; i < sx; ++i) {
      ff[i] += dx;
    }

    // convert to position correction
    for (size_t i = 0; i < sx; ++i) {
      ff[i] -= xx[i];
    }
  }
  /*
  void Reconst(const FieldCell<Scal>& uc) {
    auto sem = m.GetSem("reconst");

    // XXX: adhoc
    // overwrite u=0 if y<y0 or y>y0
    if (sem("bc-zero")) {
      Scal y0 = par->bcc_y0;
      Scal y1 = par->bcc_y1;
      auto& uu = const_cast<FieldCell<Scal>&>(uc);
      for (auto c : m.AllCells()) {
        auto x = m.GetCenter(c);
        if (x[1] < y0 || x[1] > y1) {
          uu[c] = 0.;
        }
      }
    }

    if (sem("height")) {
      fci_.Reinit(m, false);
      for (auto c : m.AllCells()) {
        Scal u = uc[c];
        if (u > 0. && u < 1.) {
          fci_[c] = true;
        }
      }
      for (auto f : m.SuFaces()) {
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        if (uc[cm] != uc[cp]) {
          fci_[cm] = true;
          fci_[cp] = true;
        }
      }
      // Compute normal and curvature [s]
      CalcNormal(uc, fci_, fc_n_, fck_);
      auto h = GetCellSize();
      // Reconstruct interface [s]
      for (auto c : m.SuCells()) {
        fc_a_[c] = R::GetLineA(fc_n_[c], uc[c], h);
      }
    }

    // Correction with normal from particles
    if (par->part && par->part_n) {
      Part(uc, sem);
      if (sem("parta")) {
        auto h = GetCellSize();
        for (auto c : m.AllCells()) {
          fc_a_[c] = R::GetLineA(fc_n_[c], uc[c], h);
        }
      }
    }
  }

*/
};
