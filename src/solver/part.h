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
  /*
  // Compute force to advance particles with exact contraints on ellipse.
  // Angle between segments linearly depends on index.
  // xx: array of positions
  // sx: size of xx
  // ll: array of lines
  // sl: size of ll
  // par->part_kattr: relaxation factor for absolute angle
  // par->part_kbend: relaxation factor for angle between segments
  // Output:
  // f: position corrections of size sx
  void PartForce2dCE(const Vect* xx, size_t sx, 
                     const std::array<Vect, 2>* ll, size_t sl, Vect* ff) {
    InterfaceForce(xx, sx, ll, sl, ff);
    Constr2(xx, sx, ff);
  }
  // Apply stretching and bending forces
  void Constr0(const Vect* xx, size_t sx, Vect* ff) {
    Vect h = GetCellSize();
    Scal hm = h.norminf();

    // stretching springs
    for (size_t i = 0; i < sx - 1; ++i) {
      Vect dx = xx[i + 1] - xx[i];
      Scal k = par->part_kstr;
      Scal d0 = par->part_h * hm; // equilibrium length
      Vect f = dx * (k * (1. - d0 / dx.norm()));
      ff[i] += f;
      ff[i + 1] -= f;
    }

    // mean angle
    Scal anm = (GetAn(xx,1) + GetAn(xx,2) + GetAn(xx,3)) / 3.; 

    // bending to mean angle
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
      Scal t = par->part_kbend * lm * lp * 
          (GetAn(xx,i) - anm * (par->part_bendmean ? 1. : 0.));
      // forces
      Vect fm = nm * (t / (2. * lm));
      Vect fp = np * (t / (2. * lp));
      // apply
      ff[im] += fm;
      ff[ip] += fp;
      ff[i] -= (fm + fp);
    }

    // scale by relaxaion factor
    for (size_t i = 0; i < sx; ++i) {
      ff[i] *= par->part_relax;
    }
  }
  // Apply exact constraints on force
  void Constr1(const Vect* xx, size_t sx, Vect* ff) {
    Vect h = GetCellSize();
    Scal hm = h.norminf();

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
    std::array<Vect, kNp> xa; // dx/dalpha
    std::array<Vect, kNp> xt; // dx/dtheta
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
    da *= par->part_kattr;
    dt *= par->part_kbend;

    // vector at angle da
    Vect ea = ra(da);
    // vector at angle dt/2
    Vect eth = ra(dt);
    // vector at angle dt
    Vect et = re(eth, eth);

    // segment vectors 
    std::array<Vect, kNp> dd;
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
    dx *= par->part_kstr;

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
  // Constraints with linear angle
  void Constr2(const Vect* xx, size_t sx, Vect* ff) {
    Vect h = GetCellSize();
    Scal hm = h.norminf();

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
    // gamma: difference for angle 1
    
    // derivatives of positions by angles
    std::array<Vect, kNp> xa; // dx/dalpha
    std::array<Vect, kNp> xt; // dx/dtheta
    std::array<Vect, kNp> xg; // dx/dgamma
    // central 
    const size_t ic = (sx - 1) / 2; 
    xa[ic] = Vect(0.);
    xt[ic] = Vect(0.);
    xg[ic] = Vect(0.);
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      Vect d;
      // forward 
      i = ic + q;
      d = xx[i] - xx[i - 1];
      xa[i] = xa[i - 1] + rr(d);
      xt[i] = xt[i - 1] + rr(d) * (q - 0.5);
      xg[i] = xg[i - 1] + rr(d) * (q - 1.);
      // backward
      i = ic - q;
      d = xx[i + 1] - xx[i];
      xa[i] = xa[i + 1] - rr(d);
      xt[i] = xt[i + 1] + rr(d) * (q - 0.5);
      xg[i] = xg[i + 1] - rr(d) * (q - 1.);
    }

    // correction of angles
    Scal da = 0.;
    Scal dt = 0.;
    Scal dg = 0.;
    for (size_t i = 0; i < sx; ++i) {
      da += ff[i].dot(xa[i]); // scale hm*hm
      dt += ff[i].dot(xt[i]);
      dg += ff[i].dot(xg[i]);
    }

    // rescale to 1
    da /= hm * hm; 
    dt /= hm * hm; 
    dg /= hm * hm; 
    
    // relaxation
    da *= par->part_kattr;
    dt *= par->part_kbend;
    dg *= par->part_kstr;

    // vector at angle da
    Vect ea = ra(da);
    // vector at angle dt/2
    Vect eth = ra(dt);
    // vector at angle dt
    Vect et = re(eth, eth);
    // vector at angle dg
    Vect eg = ra(dg);

    // segment vectors 
    std::array<Vect, kNp> dd;
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

    // apply dg
    {
      Vect e = eg;
      for (size_t q = 2; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = re(dd[i], e);
        // backward
        i = ic - q;
        dd[i] = re(dd[i], e);
        // next
        e = re(e, eg);
      }
    }

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

    // convert to position correction
    for (size_t i = 0; i < sx; ++i) {
      ff[i] -= xx[i];
    }
  }
  // Compute force to advance particles with exact constraints.
  // Assume 2d positions (x,y,0).
  // x: pointer to first particle
  // sx: size of particle string
  // l: pointer to first array of line ends
  // sl: number of lines
  // par->part_kattr: relaxation factor for absolute angle
  // par->part_kbend: relaxation factor for angle between segments
  // Output:
  // f: position corrections of size sx
  void PartForce2dC(const Vect* xx, size_t sx, 
                    const std::array<Vect, 2>* ll, size_t sl,
                    Vect* ff) {
    InterfaceForce(xx, sx, ll, sl, ff);
    Constr1(xx, sx, ff);
  }
  // Compute force to advance particles.
  // Assume 2d positions (x,y,0).
  // x: pointer to first particle
  // sx: size of particle string
  // l: pointer to first array of line ends
  // nl: number of lines
  // Output:
  // f: position corrections of size sx
  void PartForce2d(const Vect* xx, size_t sx, 
                   const std::array<Vect, 2>* ll, size_t sl, Vect* ff) {
    InterfaceForce(xx, sx, ll, sl, ff);
    Constr0(xx, sx, ff);
  }
  std::pair<Vect, Vect> GetBubble() {
    static Vect c(0);
    static Vect r(0);
    static bool ld = false;
    if (!ld) {
      std::ifstream f("../b.dat");
      f >> c[0] >> c[1] >> c[2];
      f >> r[0];
      f >> r[1];
      if (!f.good()) {
        r[1] = r[0];
      } else {
        f >> r[2];
        if (!f.good()) {
          r[2] = r[0];
        }
      }
      std::cout << "Loaded c=" << c << " r=" << r << std::endl;
      ld = true;
    }
    return std::make_pair(c, r);
  }
  void Part(const FieldCell<Scal>& uc, typename M::Sem& sem) {
    if (sem("part-comma")) {
      m.Comm(&fc_a_);
      m.Comm(&fc_n_);
    }

    bool dm = par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                            this->GetTimeStep());

    if (sem("part-seed")) {
      SeedParticles(uc);
    }

    if (sem.Nested("part-dump0")) {
      if (dm) {
        DumpParticles(0);
      }
    }

    if (sem("part-advance")) {
      // particle strings
      std::vector<IdxCell> sc; // cell
      std::vector<Vect> xx; // particle positions
      std::vector<size_t> sx; // xx index, plus element sx.size() 
      std::vector<std::array<Vect, 2>> ll; // interface lines
      std::vector<size_t> sl; // sl index, plus element sl.size() 
      std::vector<Vect> smx; // local unit x
      std::vector<Vect> smy; // local unit y
      std::vector<Vect> smc; // local center
      std::vector<Scal> sk; // curvature

      // Extract interface, project particles
      for (auto c : m.Cells()) {
        // XXX: assume fcps_[c] % kNp == 0
        for (int is = 0; is < fcps_[c] / kNp; ++is) {
          int i0 = is * kNp;
          int i1 = (is + 1) * kNp;

          sc.push_back(c);

          // Plane coordinates
          // center
          Vect rc = fcp_[c][(i0 + i1 + 1) / 2];
          // tangent, assume straight line
          Vect rt = (fcp_[c][i0 + 1] - fcp_[c][i0]);
          rt /= rt.norm();
          // interface normal, assume orthogonal to rt
          Vect n = fc_n_[c];
          n /= n.norm();
          // string plane normal
          Vect rn = rt.cross(n);
          rn /= rn.norm();

          // unit in x
          Vect mx = rt;
          // unit in y
          Vect my = n;
          // center
          Vect mc = rc;

          smx.push_back(mx);
          smy.push_back(my);
          smc.push_back(mc);

          // Transform to plane coordinates.
          // x: space coordinates
          // Returns:
          // Vect(xl, yl, 0): plane coordinates
          auto pr = [&](Vect x) -> Vect {
            Vect q(0);
            q[0] = (x - mc).dot(mx);
            q[1] = (x - mc).dot(my);
            return q;
          };

          // Copy projected particles
          sx.push_back(xx.size());
          for (int i = i0; i < i1; ++i) {
            xx.push_back(pr(fcp_[c][i]));
          }

          // Extract interface lines
          sl.push_back(ll.size());
          auto& bc = m.GetBlockCells();
          using MIdx = typename M::MIdx;
          Vect h = GetCellSize();

          const int sw = 2; // stencil halfwidth, [-sw,sw]
          const int sn = sw * 2 + 1; // stencil size

          // block of offsets
          GBlock<IdxCell, dim> bo(
              MIdx(-sw, -sw, par->dim == 2 ? 0 : -sw), 
              MIdx(sn, sn, par->dim == 2 ? 1 : sn)); 

          auto w = bc.GetMIdx(c);
          for (auto wo : bo) {
            auto cc = bc.GetIdx(w + wo);
            Scal u = uc[cc];
            Scal th = par->part_intth;
            if (fci_[c] && u > th && u < 1. - th) {
              auto xcc = m.GetCenter(cc);
              auto xx = R::GetCutPoly(xcc, fc_n_[cc], fc_a_[cc], h);
              std::array<Vect, 2> e;
              if (R::GetInterPoly(xx, rc, rn, e)) {
                // projected line ends 
                // <pncc,pe1-pe0> positively oriented
                auto pncc = pr(mc + fc_n_[cc]);
                auto pe0 = pr(e[0]);
                auto pe1 = pr(e[1]);
                if (pncc.cross_third(pe1 - pe0) < 0.) {
                  std::swap(pe0, pe1);
                }
                ll.push_back({pe0, pe1});
              }
              #if ADHOC_ATTR 
              auto cr = GetBubble();
              ll.push_back({pr(cr.first), cr.second}); 
              #endif
            }
          }
        }
      }
      sx.push_back(xx.size());
      sl.push_back(ll.size());

      assert(sx.size() == sc.size() + 1);
      assert(sl.size() == sc.size() + 1);
      assert(smx.size() == sc.size());
      assert(smy.size() == sc.size());
      assert(smc.size() == sc.size());

      sk.resize(sc.size(), 0.);

      std::vector<Vect> ff(xx.size()); // force

      // advance particles
      for (size_t it = 0; it < par->part_maxiter; ++it) {
        for (size_t i = 0; i < sc.size(); ++i) {
          if (par->part_constr == 1) {
            PartForce2dC(&(xx[sx[i]]), sx[i+1] - sx[i],
                        &(ll[sl[i]]), sl[i+1] - sl[i],
                        &(ff[sx[i]]));
          } else if (par->part_constr == 2) {
            PartForce2dCE(&(xx[sx[i]]), sx[i+1] - sx[i],
                        &(ll[sl[i]]), sl[i+1] - sl[i],
                        &(ff[sx[i]]));
          } else {
            PartForce2d(&(xx[sx[i]]), sx[i+1] - sx[i],
                        &(ll[sl[i]]), sl[i+1] - sl[i],
                        &(ff[sx[i]]));
            // freeze central particle
            Vect f = ff[(sx[i] + sx[i+1] - 1) / 2];
            for (size_t j = sx[i]; j < sx[i+1]; ++j) {
              ff[j] -= f;
            }
          }
        }

        // report error
        size_t dr = std::max<size_t>(1, 
            par->part_maxiter / par->part_report_fr);
        if (m.IsRoot() && (it % dr == 0 || it + 1 == par->part_maxiter)) {
          Vect h = GetCellSize();
          Scal hm = h.norminf();
          Scal tmax = 0.;
          Scal anmax = 0.;
          Scal anavg = 0; // average difference from mean angle
          Scal lmax = 0.; // maximum error in sement length
          size_t anavgn = 0;
          for (size_t i = 0; i < sc.size(); ++i) {
            // maximum force
            for (size_t j = sx[i]; j < sx[i+1]; ++j) {
              tmax = std::max(tmax, ff[j].norm());
            }
            Scal anm = 0.;
            size_t anmn = 0;
            // mean angle in string
            for (size_t j = sx[i] + 1; j + 1 < sx[i+1]; ++j) {
              anm += GetAn(xx.data(), j);
              ++anmn;
            }
            anm /= anmn;
            // error in angle
            for (size_t j = sx[i] + 1; j + 1 < sx[i+1]; ++j) {
              Scal e = std::abs(anm - GetAn(xx.data(), j));
              anavg += e;
              ++anavgn;
              anmax = std::max(anmax, e);
            }
            // maximum error in sement length
            for (size_t j = sx[i]; j + 1 < sx[i+1]; ++j) {
              lmax = std::max(
                  lmax, std::abs(xx[j + 1].dist(xx[j]) - par->part_h * hm));
            }
          }
          anavg /= anavgn;
          std::cout << std::setprecision(10)
              << "it=" << it 
              << " dxmax=" << tmax 
              << " anmax=" << anmax 
              << " anavg=" << anavg 
              << " lmax=" << lmax / (par->part_h * hm)
              << std::endl;
        }

        // advance
        for (size_t i = 0; i < xx.size(); ++i) {
          xx[i] += ff[i];
        }

        // copy back to field
        if (dm || it + 1 == par->part_maxiter) {
          fcps_.Reinit(m, 0);

          for (size_t i = 0; i < sc.size(); ++i) {
            IdxCell c = sc[i];
            Vect mx = smx[i];
            Vect my = smy[i];
            Vect mc = smc[i];
            for (size_t j = sx[i]; j < sx[i+1]; ++j) {
              auto q = xx[j];
              fcp_[c][fcps_[c]++] = mc + mx * q[0] + my * q[1];
            }
          }
        }

        if (dm) {
          //// XXX: disable iter dump to avoid suspender loop
          //DumpParticles(it + 1); 
        }
      }

      // compute curvature on strings
      for (size_t i = 0; i < sc.size(); ++i) {
        sk[i] = PartK(&(xx[sx[i]]), sx[i + 1] - sx[i]);
      }

      // compute curvature in cells
      fckp_.Reinit(m, 0.);
      {
        size_t i = 0;
        while (i < sc.size()) {
          IdxCell c = sc[i];
          Scal k = sk[i];
          ++i;
          size_t nk = 1;
          // average over all strings in c
          while (i < sc.size() && sc[i] == c) {
            k += sk[i];
            ++i;
            ++nk;
          }
          if (par->dim == 3) {
            k *= 2.;
          }
          fckp_[c] = k / nk;
        }
      }
      m.Comm(&fckp_);

      // compute normal
      if (par->part_n) {
        for (auto c : m.Cells()) {
          if (fcps_[c]) {
            int i = 2;
            int im = i - 1;
            int ip = i + 1;
            Vect x = fcp_[c][i];
            Vect xm = fcp_[c][im];
            Vect xp = fcp_[c][ip];
            Vect dm = x - xm;
            Vect dp = xp - x;
            Scal lm = dm.norm();
            Scal lp = dp.norm();
            // normals to segments, <nm,dm> positively oriented
            Vect nm = Vect(dm[1], -dm[0], 0.) / lm; 
            Vect np = Vect(dp[1], -dp[0], 0.) / lp;
            fc_n_[c] = (nm + np) * (-0.5);
          }
        }
      }
    }

    if (sem.Nested("part-dumplast")) {
      if (dm) {
        DumpParticles(par->part_maxiter - 1);
      }
    }
  }

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
  // Print column of datafield
  void Print(const FieldFace<Scal>& ff, std::string name) {
    using MIdx = typename M::MIdx;
    auto ibc = m.GetInBlockCells();
    auto bc = m.GetBlockCells();
    auto bf = m.GetBlockFaces();
    MIdx wb = ibc.GetBegin();
    MIdx we = ibc.GetEnd();
    MIdx wc = (wb + we - MIdx(1)) / 2;
    std::cerr << std::setw(5) << name << " = ";
    for (int i = wb[1]; i <= we[1]; ++i) {
      MIdx w = wc;
      w[1] = i;
      using Dir = typename M::Dir;
      IdxFace f = bf.GetIdx(w, Dir::j);
      std::cerr << std::setw(10) << ff[f] << " ";
    }
    std::cerr << std::endl;
  }
  // Print column of datafield
  void Print(const FieldCell<Scal>& fc, std::string name) {
    using MIdx = typename M::MIdx;
    auto ibc = m.GetInBlockCells();
    auto bc = m.GetBlockCells();
    MIdx wb = ibc.GetBegin();
    MIdx we = ibc.GetEnd();
    MIdx wc = (wb + we - MIdx(1)) / 2;
    std::cerr << std::setw(5) << name << " = ";
    for (int i = wb[1]; i < we[1]; ++i) {
      MIdx w = wc;
      w[1] = i;
      IdxCell c = bc.GetIdx(w);
      std::cerr << std::setw(10) << fc[c] << " ";
    }
    std::cerr << std::endl;
  }

  void MakeIteration() override {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      auto& uc = fc_u_.iter_curr;
      const Scal dt = this->GetTimeStep();
      for (auto c : m.Cells()) {
        uc[c] = fc_u_.time_prev[c] +  // previous time step
            dt * (*fcs_)[c]; // source
      }
    }

    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    auto& bf = m.GetBlockFaces();
    // directions, format: {dir LE, dir EI, ...}
    std::vector<size_t> dd; 
    Scal vsc; // scaling factor for ffv, used for splitting
    if (par->dim == 3) { // 3d
      if (count_ % 3 == 0) {
        dd = {0, 1, 1, 2, 2, 0};
      } else if (count_ % 3 == 1) {
        dd = {1, 2, 2, 0, 0, 1};
      } else {
        dd = {2, 0, 0, 1, 1, 2};
      }
      vsc = 0.5;
    } else {
      if (count_ % 2 == 0) {
        dd = {0, 1};
      } else {
        dd = {1, 0};
      } 
      vsc = 1.0;
    }
    for (size_t id = 0; id < dd.size(); ++id) {
      // TODO: fluxes computed twice, consider buffer
      if (sem("adv")) {
        size_t d = dd[id]; // direction as index
        Dir md(d); // direction as Dir
        MIdx wd(md); // offset in direction d
        auto& uc = fc_u_.iter_curr;
        auto& bc = m.GetBlockCells();
        auto& bf = m.GetBlockFaces();
        auto h = GetCellSize();
        auto& ffv = *ffv_; // [f]ield [f]ace [v]olume flux
        const Scal dt = this->GetTimeStep();

        ffvu_.Reinit(m);

        for (auto f : m.Faces()) {
          auto p = bf.GetMIdxDir(f);
          MIdx wf = p.first;
          Dir df = p.second;

          if (df != md) {
            continue;
          }

          // mixture flux
          const Scal v = ffv[f] * vsc;
          // upwind cell
          IdxCell cu = m.GetNeighbourCell(f, v > 0. ? 0 : 1);
          if (fci_[cu]) { // cell contains interface, flux from reconstruction
            if (id % 2 == 0) { // Euler Implicit
              // phase 2 flux
              ffvu_[f] = R::GetLineFlux(fc_n_[cu], fc_a_[cu], h, v, dt, d);
            } else { // Lagrange Explicit
              // upwind face
              IdxFace fu = bf.GetIdx(v > 0. ? wf - wd : wf + wd, md);
              // upwind mixture flux
              Scal vu = ffv[fu] * vsc;
              // phase 2 flux
              ffvu_[f] = R::GetLineFluxStr(fc_n_[cu], fc_a_[cu], h, v, vu, dt, d);
            }
          } else {
            ffvu_[f] = v * uc[cu];
          }
        }


        ffu_.Reinit(m);
        // interpolate field value to boundaries
        InterpolateB(uc, mfc_, ffu_, m);

        // override boundary upwind flux
        for (const auto& it : mfc_) {
          IdxFace f = it.GetIdx();
          CondFace* cb = it.GetValue().get(); 
          Scal v = ffv[f];
          if ((cb->GetNci() == 0) != (v > 0.)) {
            ffvu_[f] = v * ffu_[f];
            // XXX: adhoc
            // Alternating mul correction of flux
            // (done for bubble detachment)
            if (par->bcc_t0 > 0. && par->bcc_t1 > 0.) {
              Scal k0 = par->bcc_k0;
              Scal k1 = par->bcc_k1;
              Scal t0 = par->bcc_t0;
              Scal t1 = par->bcc_t1;
              Scal t = this->GetTime();
              Scal ts = t0 + t1;
              Scal ph = t / ts; 
              ph = ph - int(ph);
              ph *= ts;
              ffvu_[f] *= (ph < t0 ? k0 : k1);
            }
          }
        }

        for (auto c : m.Cells()) {
          auto w = bc.GetMIdx(c);
          const Scal lc = m.GetVolume(c);
          // faces
          IdxFace fm = bf.GetIdx(w, md);
          IdxFace fp = bf.GetIdx(w + wd, md);
          // mixture fluxes
          const Scal vm = ffv[fm] * vsc;
          const Scal vp = ffv[fp] * vsc;
          // mixture cfl
          const Scal sm = vm * dt / lc;
          const Scal sp = vp * dt / lc;
          const Scal ds = sp - sm;
          // phase 2 fluxes
          Scal qm = ffvu_[fm];
          Scal qp = ffvu_[fp];
          // phase 2 cfl
          const Scal lm = qm * dt / lc;
          const Scal lp = qp * dt / lc;
          const Scal dl = lp - lm;
          if (id % 2 == 0) { // Euler Implicit
            uc[c] = (uc[c] - dl) / (1. - ds);
          } else { // Lagrange Explicit
            uc[c] = uc[c] * (1. + ds) - dl;
          }
        }

        // clip
        const Scal th = par->clipth;
        for (auto c : m.Cells()) {
          Scal& u = uc[c];
          if (u < th) {
            u = 0.;
          } else if (u > 1. - th) {
            u = 1.;
          }
        }
        m.Comm(&uc);
      }
      if (sem.Nested("reconst")) {
        Reconst(fc_u_.iter_curr);
      }
    }

    // Curvature from gradient of volume fraction
    if (par->curvgrad && sem("curv")) {
      ffu_ = Interpolate(fc_u_.iter_curr, mfc_, m); // [s]
      fcg_ = Gradient(ffu_, m); // [s]
      ffg_ = Interpolate(fcg_, mfvz_, m); // [i]

      fck_.Reinit(m, GetNan<Scal>()); // curvature [i]
      for (auto c : m.Cells()) {
        if (!fci_[c]) {
          continue;
        }
        Scal s = 0.;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          auto& g = ffg_[f];
          auto n = g / g.norm();  // inner normal
          s += -n.dot(m.GetOutwardSurface(c, q));
        }
        fck_[c] = s / m.GetVolume(c);
      }
    }

    if (sem("curvcomm")) {
      m.Comm(&fck_);
    }

    if (par->dumppoly) {
      bool dm = par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                              this->GetTimeStep());
      if (dm && sem("dumppoly")) {
        DumpPoly();
      }
    }

    if (par->part) {
      Part(fc_u_.iter_curr, sem);
    }

    if (sem("stat")) {
      this->IncIter();
      ++count_;
    }
  }
  void FinishStep() override {
    fc_u_.time_curr = fc_u_.iter_curr;
    this->IncTime();
  }
  const FieldCell<Scal>& GetField(Layers l) const override {
    return fc_u_.Get(l);
  }
  const FieldCell<Scal>& GetAlpha() const {
    return fc_a_;
  }
  const FieldCell<Vect>& GetNormal() const {
    return fc_n_;
  }
  const FieldCell<Scal>& GetCurv() const override {
    return par->part_k ? fckp_ : fck_;
  }
  // curvature from height function
  const FieldCell<Scal>& GetCurvH() const {
    return fck_;
  }
  // curvature from particles
  const FieldCell<Scal>& GetCurvP() const {
    return fckp_;
  }
  using P::GetField;
};

} // namespace solver
*/
};
