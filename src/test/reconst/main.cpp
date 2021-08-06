// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "func/init_u.h"
#include "overlap/overlap.h"
#include "solver/vof.h"

using Scal = double;
using Vect = generic::Vect<Scal, 3>;
using Vect2 = generic::Vect<Scal, 2>;
using R = Reconst<Scal>;

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

void TestRandom() {
  std::cerr << "Check u=u(a(u)) on random input" << std::endl;
  auto t = [](Vect msk, std::string s) {
    std::default_random_engine g(0);
    std::uniform_real_distribution<double> f(0., 1.);
    size_t ni = 10000;
    Scal e = 0.; // error
    Scal em = 0.; // error max
    size_t nis = 0; // sucessful samples
    for (size_t i = 0; i < ni; ++i) {
      Vect n(f(g), f(g), f(g));
      n *= msk;
      n /= n.norm();
      Scal u = f(g);

      Scal a = R::GetLineA1(n, u);
      Scal uu = R::GetLineU1(n, a);
      e += std::abs(u - uu);
      em = std::max(em, std::abs(u - uu));

      /*
      std::cout << "n=" << n
        << " u=" << u
        << " uu=" << uu
        << " a=" << a
        << std::endl;
      */

      /*

ba:
set n=(1,1,0)
n=(0.7071067811865475,0.7071067811865475,0) u=0.5 h=(1,1,1) dx=1
a=-0.2071067811865476 dv=0.2499999999999999 dve=0.5

ch:
set n=(1,1,0)
n=(0.7071067811865475,0.7071067811865475,0) u=0.5 h=(1,1,1) dx=1
a=9.93410742555767e-09 dv=0.5000000140489494 dve=0.5

       * */

      Vect h(0.1, 0.2, 0.3);
      e += std::abs(R::GetLineU(n, R::GetLineA(n, u, h), h) - u);

      ++nis;
    }
    e /= nis;
    std::cerr << s + ": eavg=" << e << " emax=" << em << " samples=" << ni
              << std::endl;
    // assert(e < 1e-14);
  };
  t(Vect(1., 0., 0.), "x");
  t(Vect(0., 1., 0.), "y");
  t(Vect(0., 0., 1.), "z");
  t(Vect(1., 1., 0.), "xy");
  t(Vect(0., 1., 1.), "yz");
  t(Vect(1., 0., 1.), "zx");
  t(Vect(1., 1., 1.), "xyz");
}

// rectangular cell test
void TestRect() {
  auto p = [](Scal nx, Scal ny, Scal nz, Scal u, Scal hx, Scal hy, Scal hz) {
    Vect n(nx, ny, nz);
    Vect h(hx, hy, hz);
    n /= n.norm();
    Scal a = R::GetLineA(n, u, h);
    Scal ua = R::GetLineU(n, a, h);
    Scal dx = 0.5 * hx; // advection in x
    Vect hs(dx, hy, 1.); // size of acceptor
    Vect d = Vect(hx * 0.5 + dx * 0.5, 0., 0.); // shift of center
    Scal uadv = R::GetLineU(n, a - n.dot(d), hs);
    std::cerr << "n=" << n << " u=" << u << " h=" << h << " a=" << a
              << " u(a)=" << ua << " voladv=" << uadv * hs.prod() << std::endl;
  };
  auto v = [](Scal nx, Scal ny, Scal nz, Scal a, Scal ue) {
    Vect n(nx, ny, nz);
    Scal u = R::GetLineU(n, a, Vect(1.));
    std::cerr << "n=" << n << " a=" << a << " u=" << u << " ue=" << ue
              << std::endl;
    assert(std::abs(u - ue) < 1e-12);
  };
  std::cerr << "check rectangular cell" << std::endl;
  p(1., 1., 0., 0., 0.5, 1., 1.);
  p(1., 1., 0., 0.5, 0.1, 1., 1.);
  p(1., 1., 0., 0.25, 1., 1., 1.);
  p(1., 1., 0., 0.25, 0.1, 1., 1.);
  p(0., 1., 0., 0.5, 1., 1., 1.);
  v(0., 1., 0., 0., 0.5);
  v(1., 1., 1., -Vect(0.5).norm1(), 0.);
  v(1., 1., 1., Vect(0.5).norm1(), 1.);
  v(1., 1., 1., -0.5, 1. / 6.);
}

// volume surplus test
void TestVol() {
  // nx, ny: normal
  // u: volume fraction
  // dx: displacement
  // dve: volume surplus (exact)
  // xx: direction x
  auto fg = [](Scal nx, Scal ny, Scal u, Scal hx, Scal hy, Scal dx, Scal dve,
               bool xx) {
    Vect n(nx, ny, 0.);
    Vect h(hx, hy, 1.);
    n /= n.norm();
    Scal a = R::GetLineA(n, u, h);
    Scal dv;
    dv = R::GetLineVol(n, a, h, dx, xx ? 0 : 1);
    std::cerr << std::setprecision(16) << "n=" << n << " u=" << u << " h=" << h
              << " dx=" << dx << " a=" << a << " dv=" << dv << " dve=" << dve
              << std::endl;
    assert(std::abs(dv - dve) < 1e-6);
  };
  auto fy = [&fg](
                Scal nx, Scal ny, Scal u, Scal hx, Scal hy, Scal dx, Scal dve) {
    fg(nx, ny, u, hx, hy, dx, dve, false);
  };
  // f=fx
  auto f = [&fg](
               Scal nx, Scal ny, Scal u, Scal hx, Scal hy, Scal dx, Scal dve) {
    fg(nx, ny, u, hx, hy, dx, dve, true);
  };

  Scal nx, ny, u, hx, hy, dx;
  Scal nx0, ny0, u0, hx0, hy0, dx0;

  std::cerr << "check volume surplus" << std::endl;
  nx0 = 1.;
  ny0 = 0.;
  u0 = 0.5;
  hx0 = 1.;
  hy0 = 1.;
  dx0 = 0.;
  nx = nx0;
  ny = ny0;
  u = u0;
  hx = hx0;
  hy = hy0;
  dx = dx0;

  std::cerr << "base" << std::endl;
  f(nx, ny, u, hx, hy, dx = 1., 0.5);
  f(nx, ny, u, hx, hy, dx = 0.5, 0.);
  f(nx, ny, u, hx, hy, dx = -1., 0.5);
  f(nx, ny, u, hx, hy, dx = -0.5, 0.5);

  std::cerr << "set u=1" << std::endl;
  u = 1.;
  f(nx, ny, u, hx, hy, dx = 0., 0.);
  f(nx, ny, u, hx, hy, dx = 1., 1.);
  f(nx, ny, u, hx, hy, dx = 0.5, 0.5);
  u = u0;

  std::cerr << "set n=(-1,0,0)" << std::endl;
  nx = -1.;
  f(nx, ny, u, hx, hy, dx = 1., 0.5);
  f(nx, ny, u, hx, hy, dx = 0.5, 0.5);
  f(nx, ny, u, hx, hy, dx = -1., 0.5);
  f(nx, ny, u, hx, hy, dx = -0.5, 0.);
  nx = nx0;

  std::cerr << "set hx=0.5" << std::endl;
  hx = 0.5;
  f(nx, ny, u, hx, hy, dx = 0.5, 0.25);
  f(nx, ny, u, hx, hy, dx = 0.25, 0.);
  f(nx, ny, u, hx, hy, dx = -0.5, 0.25);
  f(nx, ny, u, hx, hy, dx = -0.25, 0.25);
  hx = hx0;

  std::cerr << "set n=(0,1,0)" << std::endl;
  nx = 0.;
  ny = 1.;
  f(nx, ny, u, hx, hy, dx = 1., 0.5);
  f(nx, ny, u, hx, hy, dx = 0.5, 0.25);
  f(nx, ny, u, hx, hy, dx = -1., 0.5);
  f(nx, ny, u, hx, hy, dx = -0.5, 0.25);
  nx = nx0;
  ny = ny0;

  std::cerr << "set n=(1,1,0)" << std::endl;
  nx = 1.;
  ny = 1.;
  f(nx, ny, u, hx, hy, dx = 1., 0.5);
  f(nx, ny, u, hx, hy, dx = 0.5, 0.125);
  f(nx, ny, u, hx, hy, dx = -1., 0.5);
  f(nx, ny, u, hx, hy, dx = -0.5, 0.375);
  nx = nx0;
  ny = ny0;

  std::cerr << "set u=0" << std::endl;
  u = 0.;
  f(nx, ny, u, hx, hy, dx = 1., 0.);
  f(nx, ny, u, hx, hy, dx = 0.5, 0.);
  f(nx, ny, u, hx, hy, dx = -1., 0.);
  f(nx, ny, u, hx, hy, dx = -0.5, 0.);
  nx = -nx;
  f(nx, ny, u, hx, hy, dx = 1., 0.);
  f(nx, ny, u, hx, hy, dx = 0.5, 0.);
  f(nx, ny, u, hx, hy, dx = -1., 0.);
  f(nx, ny, u, hx, hy, dx = -0.5, 0.);
  u = u0;
  nx = nx0;

  std::cerr << "set dx=0.01" << std::endl;
  hx = 0.1;
  hy = 0.1;
  f(nx, ny, u = 0.6, hx, hy, dx = 0.01, 0.);
  f(nx, ny, u = 0.55, hx, hy, dx = 0.01, 0.);
  f(nx, ny, u = 0.4, hx, hy, dx = 0.01, 0.);
  f(ny, nx, u = 0.6, hx, hy, dx = 0.01, 0.006 * hy);
  f(ny, nx, u = 0.4, hx, hy, dx = 0.01, 0.004 * hy);
  fy(nx, ny, u = 0.6, hx, hy, dx = 0.01, 0.006 * hx);
  fy(nx, ny, u = 0.4, hx, hy, dx = 0.01, 0.004 * hx);
  fy(ny, nx, u = 0.6, hx, hy, dx = 0.01, 0.);
  fy(ny, nx, u = 0.4, hx, hy, dx = 0.01, 0.);
  u = u0;
  nx = nx0;
  hx = hx0;
  hy = hy0;
}

// stretching volume surplus test
void TestVolStr() {
  // nx, ny: normal
  // u: volume fraction
  // dx: displacement
  // dxu: displacement upwind
  // dve: volume surplus (exact)
  auto f = [](Scal nx, Scal ny, Scal u, Scal hx, Scal hy, Scal dx, Scal dxu,
              Scal dve) {
    Vect n(nx, ny, 0.);
    Vect h(hx, hy, 1.);
    n /= n.norm();
    Scal a = R::GetLineA(n, u, h);
    Scal dv;
    dv = R::GetLineVolStr(n, a, h, dx, dxu, 0);
    std::cerr << "n=" << n << " u=" << u << " h=" << h << " dx=" << dx
              << " dxu=" << dxu << " a=" << a << " dv=" << dv << " dve=" << dve
              << std::endl;
    assert(std::abs(dv - dve) < 1e-7);
  };

  Scal nx, ny, u, hx, hy, dx, dxu;
  Scal nx0, ny0, u0, hx0, hy0, dx0, dxu0;

  std::cerr << "check volume surplus" << std::endl;
  nx0 = 1.;
  ny0 = 0.;
  u0 = 0.5;
  hx0 = 1.;
  hy0 = 1.;
  dx0 = 0.;
  dxu0 = 0.;
  nx = nx0;
  ny = ny0;
  u = u0;
  hx = hx0;
  hy = hy0;
  dx = dx0;
  dxu = dxu0;

  std::cerr << "base" << std::endl;
  f(nx, ny, u, hx, hy, dx = 1., dxu = 1., 0.5);
  f(nx, ny, u, hx, hy, dx = -1., dxu = -1., 0.5);
  f(nx, ny, u, hx, hy, dx = 1., dxu = 0., 0.);
  f(nx, ny, u, hx, hy, dx = -1., dxu = 0., 1.);
  f(nx, ny, u, hx, hy, dx = 0.5, dxu = 0., 0.);
  f(nx, ny, u, hx, hy, dx = -1.5, dxu = -1., 0.75);
  f(nx, ny, u, hx, hy, dx = -0.5, dxu = 0., 0.5);

  std::cerr << "set n=(0,1,0)" << std::endl;
  nx = 0.;
  ny = 1.;
  f(nx, ny, u, hx, hy, dx = 1., dxu = 1., 0.5);
  f(nx, ny, u, hx, hy, dx = -1., dxu = -1., 0.5);
  f(nx, ny, u, hx, hy, dx = 1., dxu = 0., 0.5);
  f(nx, ny, u, hx, hy, dx = -1., dxu = 0., 0.5);
  f(nx, ny, u, hx, hy, dx = 0.5, dxu = 0.5, 0.25);
  f(nx, ny, u, hx, hy, dx = -0.5, dxu = -0.5, 0.25);
  nx = nx0;
  ny = ny0;

  std::cerr << "set n=(1,1,0)" << std::endl;
  nx = 1.;
  ny = 1.;
  f(nx, ny, u, hx, hy, dx = 1., dxu = 1., 0.5);
  f(nx, ny, u, hx, hy, dx = -1., dxu = -1., 0.5);
  f(nx, ny, u, hx, hy, dx = 1., dxu = 0., 0.25);
  f(nx, ny, u, hx, hy, dx = -1., dxu = 0., 0.75);
  nx = nx0;
  ny = ny0;
}

void Plot() {
  Vect n(0.1, 0.2, 0.29);
  {
    std::string fn = "vof_a.dat";
    std::ofstream f(fn);
    for (Scal u = 0.; u < 1. + 1e-6; u += 0.001) {
      f << u << " " << R::GetLineA1(n, u) << "\n";
    }
  }

  {
    std::string fn = "vof_u.dat";
    std::ofstream f(fn);
    Scal nh = n.norm1() * 0.5;
    for (Scal a = -nh; a < nh + 1e-6; a += 0.001 * nh * 2.) {
      f << a << " " << R::GetLineU1(n, a) << "\n";
    }
  }
}

void TestFit() {
  std::cerr << "check plane fitting" << std::endl;

  using V = Vect;
  using VV = generic::Vect<V, 3>;

  auto f = [](const VV& a, const V& xe) {
    auto x = R::SolveSingular(a);
    std::cout << "a=" << a << " x=" << x << " xe=" << xe << std::endl;
    assert(x.dist(xe) < 1e-12);
  };

  f(VV(V(-2., 1., 1.), V(1., -2., 1.), V(1., 1., -2.)), V(1., 1., 1.));
  f(VV(V(1., 0., 0.), V(0., 1., 0.), V(1., 1., 0.)), V(0., 0., 1.));

  auto g = [](const std::vector<V>& xx, const V& ne) {
    auto n = R::GetFitN(xx);
    std::cout << "xx=" << xx << " n=" << n << " ne=" << ne << std::endl;
    assert(n.dist(ne) < 1e-12);
  };

  g({V(0., 0., 0.), V(0., 1., 0.), V(0., 1., 1.), V(0., 0., 1.)},
    V(1., 0., 0.));
}

void TestNearest() {
  std::cerr << "check nearest to polygon" << std::endl;
  using V = Vect;
  auto f = [](V x, const std::vector<V>& xx, V n, V xne) {
    auto xn = R::GetNearest(x, xx, n);
    std::cout << "tri:"
              << " x=" << x << " xx=" << xx << " n=" << n << " xn=" << xn
              << " xne=" << xne << std::endl;
    assert(xne.dist(xn) < 1e-12);
  };

  std::vector<V> xx = {V(0., 0., 0.), V(2., 0., 0.), V(0., 2., 0.)};
  V n(0., 0., 6.);
  f(V(10., 10., 1.), xx, n, V(1., 1., 0.));
  f(V(10., 0., 1.), xx, n, V(2., 0., 0.));
  f(V(0., 10., 1.), xx, n, V(0., 2., 0.));
  f(V(0.1, 0.2, 1.), xx, n, V(0.1, 0.2, 0.));
}

// test init from level-set
void TestLevelSet() {
  std::cerr << "check volume from level-set" << std::endl;
  auto f = [](const Vect& x) -> Scal { return 1. - x.dot(x); };

  auto t = [&](const Vect& xc, const Vect& h, Scal ue) {
    Scal u = GetLevelSetVolume<Scal, 3>(f, xc, h);
    std::cout << "circle:"
              << " xc=" << xc << " h=" << h << " u=" << u << " ue=" << ue
              << std::endl;
    assert(std::abs(u - ue) < 1e-12);
  };

  t(Vect(1., 0., 0.), Vect(0.1), 0.5);
  t(Vect(1.1, 0., 0.), Vect(0.1), 0.);
  t(Vect(0.9, 0., 0.), Vect(0.1), 1.);
  t(Vect(1. / sqrt(2.), 1. / sqrt(2.), 0.), Vect(0.1), 0.5);
}

// intersection of line and convex plane polygon
void TestInter() {
  std::cerr << "check Reconst::GetInter" << std::endl;
  using V = Vect2;
  // xx: polygon
  // xc: point on line
  // t: line tangent
  // aae: reference intersection
  auto f = [](const std::vector<V>& xx, V xc, V t, std::array<Scal, 2> aae) {
    std::array<Scal, 2> aa;
    R::GetInter(xx.data(), xx.size(), xc, t, aa);
    if (aa[1] < aa[0]) {
      std::swap(aa[0], aa[1]);
    }

    std::cout << "GetInter:"
              << " xx=" << xx << " xc=" << xc << " t=" << t << " aa=(" << aa[0]
              << "," << aa[1] << ")"
              << " aae=(" << aae[0] << "," << aae[1] << ")" << std::endl;
    assert(std::abs(aa[0] - aae[0]) < 1e-12);
    assert(std::abs(aa[1] - aae[1]) < 1e-12);
  };

  std::vector<V> xx = {V(0., 0.), V(1., 0.), V(0., 1.)};
  f(xx, V(0., 0.), V(1., 1.), {0., 0.5});
  f(xx, V(-1, -1.), V(1., 1.), {1., 1.5});
  f(xx, V(0., 0.5), V(1., 0.), {0., 0.5});
}

void TestOverlap() {
  Vect h(0.1);
  Vect c(0.);
  Scal r = 1.;
  std::cout << GetSphereOverlap(Vect(0), h, c, r) << std::endl;
  std::cout << GetSphereOverlap(Vect(1), h, c, r) << std::endl;
  std::cout << GetSphereOverlap(Vect(1., 0., 0.), h, c, r) << std::endl;
  std::cout << GetSphereOverlap(Vect(0.7, 0.7, 0.), h, c, r) << std::endl;
}

void TestArea() {
  std::cerr << "check area of polygon" << std::endl;
  using V = Vect;

  // xx: polygon
  // ae: area exact
  auto f = [](const std::vector<V>& xx, Scal ae) {
    auto xc = R::GetCenter(xx);
    auto n = (xx[0] - xc).cross(xx[1] - xc);
    Scal a = R::GetArea(xx, n);
    std::cout << "area:"
              << " xx=" << xx << " n=" << n << " a=" << a << " ae=" << ae
              << std::endl;
    assert(std::abs(a - ae) < 1e-12);
  };

  {
    f({V(0., 0., 0.), V(1., 0., 0.), V(0., 1., 0.)}, 0.5);
    f({V(0., 0., 0.), V(0., 1., 0.), V(0., 0., 1.)}, 0.5);
    f({V(0., 5., 0.), V(0., 6., 0.), V(0., 5., 1.)}, 0.5);
    f({V(0., 0., 0.), V(0., 1., 10.), V(0., 0., 1.)}, 0.5);
    f({V(0., 0., 0.), V(1., 1., 0.), V(0., 0., 1.)}, std::sqrt(2.) * 0.5);
    f({V(0., 0., 0.), V(1., 1., 0.), V(1., 1., 0.), V(0., 0., 1.)},
      std::sqrt(2.) * 0.5);
    f({V(0., 0., 0.), V(1., 0., 0.), V(1., 1., 0.), V(0., 1., 0.)}, 1.);
    // Rotation around axis d at angle p
    auto r = [](Vect x, size_t d, Scal p) {
      size_t u = (d + 1) % 3;
      size_t v = (d + 2) % 3;
      Scal c = std::cos(p);
      Scal s = std::sin(p);
      Vect xr = x;
      xr[u] = c * x[u] - s * x[v];
      xr[v] = s * x[u] + c * x[v];
      return xr;
    };
    auto rv = [&r](std::vector<Vect> xx, size_t d, Scal p) {
      for (auto& x : xx) {
        x = r(x, d, p);
      }
      return xx;
    };
    // Scale by k
    auto sv = [](std::vector<Vect> xx, Scal k) {
      for (auto& x : xx) {
        x *= k;
      }
      return xx;
    };
    auto xx = {V(0., 0., 0.), V(1., 0., 0.), V(1., 1., 0.), V(0., 1., 0.)};
    f(rv(rv(rv(xx, 0, 0.5), 1, 0.9), 2, 4.3), 1.);
    f(sv(rv(rv(rv(xx, 0, 0.5), 1, 0.9), 2, 4.3), 7.), 49.);
  }
}

int main() {
  TestArea();
  TestOverlap();
  TestInter();
  TestLevelSet();
  TestNearest();
  TestFit();
  Plot();
  TestRect();
  TestRandom();
  TestVol();
  TestVolStr();
}
