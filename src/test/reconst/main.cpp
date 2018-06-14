#undef NDEBUG
#include <iostream>
#include <string>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <functional>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <random>

#include "solver/vof.h"

using Scal = double;
using Vect = GVect<Scal, 3>;

using solver::cube;

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
    std::uniform_real_distribution<double> f(0.,1.);
    size_t ni = 10000;
    Scal e = 0.; // error
    Scal em = 0.; // error max
    size_t nis = 0; // sucessful samples
    for (size_t i = 0; i < ni; ++i) {
      Vect n(f(g), f(g), f(g));
      n *= msk;
      n /= n.norm();
      Scal u = f(g);
      using namespace solver;

      Scal a = GetLineA1(n, u);
      Scal uu = GetLineU1(n, a);
      e += std::abs(u - uu);
      em = std::max(em, std::abs(u - uu));

      /*
      std::cout << "n=" << n 
        << " u=" << u 
        << " uu=" << uu 
        << " a=" << a 
        << std::endl;
      */

      Vect h(0.1, 0.2, 0.3);
      e += std::abs(solver::GetLineU(n, solver::GetLineA(n, u, h), h) - u);

      ++nis;
    }
    e /= nis;
    std::cerr << s + ": eavg=" << e 
        << " emax=" << em
        << " samples=" << ni 
        << std::endl;
    //assert(e < 1e-14);
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
    Scal a = solver::GetLineA(n, u, h);
    Scal ua = solver::GetLineU(n, a, h);
    Scal dx = 0.5 * hx; // advection in x
    Vect hs(dx, hy, 1.); // size of acceptor
    Vect d = Vect(hx * 0.5 + dx * 0.5, 0., 0.); // shift of center
    Scal uadv = solver::GetLineU(n, a - n.dot(d), hs); 
    std::cerr 
        << "n=" << n
        << " u=" << u
        << " h=" << h
        << " a=" << a
        << " u(a)=" << ua
        << " voladv=" << uadv * hs.prod()
        << std::endl;
  };
  auto v = [](Scal nx, Scal ny, Scal nz, Scal a, Scal ue) {
    Vect n(nx, ny, nz);
    Scal u = solver::GetLineU(n, a, Vect(1.));
    std::cerr 
        << "n=" << n
        << " a=" << a
        << " u=" << u
        << " ue=" << ue
        << std::endl;
    assert(std::abs(u - ue) < 1e-12);
  };
  std::cerr << "check rectangular cell" << std::endl;
  p(1., 1., 0., 0.,  0.5 , 1., 1.);
  p(1., 1., 0., 0.5 , 0.1, 1., 1.);
  p(1., 1., 0., 0.25 , 1., 1., 1.);
  p(1., 1., 0., 0.25 , 0.1, 1., 1.);
  p(0., 1., 0., 0.5 , 1., 1., 1.);
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
  auto f = [](Scal nx, Scal ny, Scal u, Scal hx, Scal hy, Scal dx, Scal dve, bool xx=true) {
    Vect n(nx, ny, 0.);
    Vect h(hx, hy, 1.);
    n /= n.norm();
    Scal a = solver::GetLineA(n, u, h);
    Scal dv;
    dv = solver::GetLineVol(n, a, h, dx, xx ? 0 : 1);
    std::cerr << std::setprecision(16)
        << "n=" << n
        << " u=" << u
        << " h=" << h
        << " dx=" << dx
        << " a=" << a
        << " dv=" << dv
        << " dve=" << dve
        << std::endl;
    assert(std::abs(dv - dve) < 1e-6);
  };

  Scal nx, ny, u, hx, hy, dx;
  Scal nx0, ny0, u0, hx0, hy0, dx0;

  std::cerr << "check volume surplus" << std::endl;
  nx0 = 1.; ny0 = 0.; u0 = 0.5; hx0 = 1.; hy0 = 1.; dx0 = 0.;
  nx = nx0; ny = ny0; u = u0; hx = hx0; hy = hy0; dx = dx0;

  std::cerr << "base" << std::endl;
  f(nx, ny, u, hx, hy, dx=1., 0.5);
  f(nx, ny, u, hx, hy, dx=0.5, 0.);
  f(nx, ny, u, hx, hy, dx=-1., 0.5);
  f(nx, ny, u, hx, hy, dx=-0.5, 0.5);

  std::cerr << "set u=1" << std::endl;
  u = 1.;
  f(nx, ny, u, hx, hy, dx=0., 0.);
  f(nx, ny, u, hx, hy, dx=1., 1.);
  f(nx, ny, u, hx, hy, dx=0.5, 0.5);
  u = u0;

  std::cerr << "set n=(-1,0,0)" << std::endl;
  nx = -1.;
  f(nx, ny, u, hx, hy, dx=1., 0.5);
  f(nx, ny, u, hx, hy, dx=0.5, 0.5);
  f(nx, ny, u, hx, hy, dx=-1., 0.5);
  f(nx, ny, u, hx, hy, dx=-0.5, 0.);
  nx = nx0;

  std::cerr << "set hx=0.5" << std::endl;
  hx = 0.5; 
  f(nx, ny, u, hx, hy, dx=0.5, 0.25);
  f(nx, ny, u, hx, hy, dx=0.25, 0.);
  f(nx, ny, u, hx, hy, dx=-0.5, 0.25);
  f(nx, ny, u, hx, hy, dx=-0.25, 0.25);
  hx = hx0;

  std::cerr << "set n=(0,1,0)" << std::endl;
  nx = 0.; ny = 1.;
  f(nx, ny, u, hx, hy, dx=1., 0.5);
  f(nx, ny, u, hx, hy, dx=0.5, 0.25);
  f(nx, ny, u, hx, hy, dx=-1., 0.5);
  f(nx, ny, u, hx, hy, dx=-0.5, 0.25);
  nx = nx0; ny = ny0;

  std::cerr << "set n=(1,1,0)" << std::endl;
  nx = 1.; ny = 1.;
  f(nx, ny, u, hx, hy, dx=1., 0.5);
  f(nx, ny, u, hx, hy, dx=0.5, 0.125);
  f(nx, ny, u, hx, hy, dx=-1., 0.5);
  f(nx, ny, u, hx, hy, dx=-0.5, 0.375);
  nx = nx0; ny = ny0;

  std::cerr << "set u=0" << std::endl;
  u = 0.;
  f(nx, ny, u, hx, hy, dx=1., 0.);
  f(nx, ny, u, hx, hy, dx=0.5, 0.);
  f(nx, ny, u, hx, hy, dx=-1., 0.);
  f(nx, ny, u, hx, hy, dx=-0.5, 0.);
  nx = -nx;
  f(nx, ny, u, hx, hy, dx=1., 0.);
  f(nx, ny, u, hx, hy, dx=0.5, 0.);
  f(nx, ny, u, hx, hy, dx=-1., 0.);
  f(nx, ny, u, hx, hy, dx=-0.5, 0.);
  u = u0;
  nx = nx0;

  std::cerr << "set dx=0.01" << std::endl;
  hx = 0.1;
  hy = 0.1;
  f(nx, ny, u=0.6, hx, hy, dx=0.01, 0.);
  f(nx, ny, u=0.55, hx, hy, dx=0.01, 0.);
  f(nx, ny, u=0.4, hx, hy, dx=0.01, 0.);
  f(ny, nx, u=0.6, hx, hy, dx=0.01, 0.006 * hy);
  f(ny, nx, u=0.4, hx, hy, dx=0.01, 0.004 * hy);
  f(nx, ny, u=0.6, hx, hy, dx=0.01, 0.006 * hx, false);
  f(nx, ny, u=0.4, hx, hy, dx=0.01, 0.004 * hx, false);
  f(ny, nx, u=0.6, hx, hy, dx=0.01, 0., false);
  f(ny, nx, u=0.4, hx, hy, dx=0.01, 0., false);
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
  auto f = [](Scal nx, Scal ny, Scal u, Scal hx, Scal hy, 
              Scal dx, Scal dxu, Scal dve, bool xx=true) {
    Vect n(nx, ny, 0.);
    Vect h(hx, hy, 1.);
    n /= n.norm();
    Scal a = solver::GetLineA(n, u, h);
    Scal dv;
    dv = solver::GetLineVolStr(n, a, h, dx, dxu, xx ? 0 : 1);
    std::cerr 
        << "n=" << n
        << " u=" << u
        << " h=" << h
        << " dx=" << dx
        << " dxu=" << dxu
        << " a=" << a
        << " dv=" << dv
        << " dve=" << dve
        << std::endl;
    assert(std::abs(dv - dve) < 1e-7);
  };

  Scal nx, ny, u, hx, hy, dx, dxu;
  Scal nx0, ny0, u0, hx0, hy0, dx0, dxu0;

  std::cerr << "check volume surplus" << std::endl;
  nx0 = 1.; ny0 = 0.; u0 = 0.5; hx0 = 1.; hy0 = 1.; dx0 = 0.; dxu0 = 0.;
  nx = nx0; ny = ny0; u = u0; hx = hx0; hy = hy0; dx = dx0; dxu = dxu0;

  std::cerr << "base" << std::endl;
  f(nx, ny, u, hx, hy, dx=1., dxu=1., 0.5);
  f(nx, ny, u, hx, hy, dx=-1., dxu=-1., 0.5);
  f(nx, ny, u, hx, hy, dx=1., dxu=0., 0.);
  f(nx, ny, u, hx, hy, dx=-1., dxu=0., 1.);
  f(nx, ny, u, hx, hy, dx=0.5, dxu=0., 0.);
  f(nx, ny, u, hx, hy, dx=-1.5, dxu=-1., 0.75);
  f(nx, ny, u, hx, hy, dx=-0.5, dxu=0., 0.5);

  std::cerr << "set n=(0,1,0)" << std::endl;
  nx = 0.; ny = 1.;
  f(nx, ny, u, hx, hy, dx=1., dxu=1., 0.5);
  f(nx, ny, u, hx, hy, dx=-1., dxu=-1., 0.5);
  f(nx, ny, u, hx, hy, dx=1., dxu=0., 0.5);
  f(nx, ny, u, hx, hy, dx=-1., dxu=0., 0.5);
  f(nx, ny, u, hx, hy, dx=0.5, dxu=0.5, 0.25);
  f(nx, ny, u, hx, hy, dx=-0.5, dxu=-0.5, 0.25);
  nx = nx0; ny = ny0;

  std::cerr << "set n=(1,1,0)" << std::endl;
  nx = 1.; ny = 1.;
  f(nx, ny, u, hx, hy, dx=1., dxu=1., 0.5);
  f(nx, ny, u, hx, hy, dx=-1., dxu=-1., 0.5);
  f(nx, ny, u, hx, hy, dx=1., dxu=0., 0.25);
  f(nx, ny, u, hx, hy, dx=-1., dxu=0., 0.75);
  nx = nx0; ny = ny0;
}

// center of line test
void TestCenter() {
  // nx, ny: normal
  // u: volume fraction
  // hx, hy: cell size
  // ce: center exact
  auto f = [](Scal nx, Scal ny, Scal u, Scal hx, Scal hy, Vect ce) {
    Vect n(nx, ny, 0.);
    Vect h(hx, hy, 1.);
    n /= n.norm();
    Scal a = solver::GetLineA(n, u, h);
    std::array<Vect, 2> xx = solver::GetLineEnds(n, a, h); 
    Vect c = solver::GetLineC(n, a, h);
    std::cerr 
        << "n=" << n
        << " u=" << u
        << " h=" << h
        << " a=" << a
        << " x0=" << xx[0]
        << " x1=" << xx[1]
        << " c=" << c
        << " ce=" << ce
        << std::endl;
    assert((c - ce).norm() < 1e-8);
  };

  Scal nx, ny, u, hx, hy;
  Scal nx0, ny0, u0, hx0, hy0;

  std::cerr << "check line center" << std::endl;
  nx0 = 1.; ny0 = 0.; u0 = 0.5; hx0 = 1.; hy0 = 1.; 
  nx = nx0; ny = ny0; u = u0; hx = hx0; hy = hy0;

  std::cerr << "base" << std::endl;
  f(nx, ny, u, hx, hy, Vect(0., 0., 0.));
  f(ny, nx, u, hx, hy, Vect(0., 0., 0.));
  f(nx, ny, u=0.25, hx, hy, Vect(-0.25, 0., 0.));
  f(ny, nx, u=0.25, hx, hy, Vect(0., -0.25, 0.));
  f(nx, nx, u=u0, hx, hy, Vect(0., 0., 0.));
}

// nearest point to line
void TestLineNearest() {
  // x: target point
  // nx, ny: normal
  // u: volume fraction
  // hx, hy: cell size
  // ce: center exact
  auto f = [](Vect x, Scal nx, Scal ny, Scal u, Scal hx, Scal hy, Vect pe) {
    Vect n(nx, ny, 0.);
    Vect h(hx, hy, 1.);
    n /= n.norm();
    Scal a = solver::GetLineA(n, u, h);
    std::array<Vect, 2> xx = solver::GetLineEnds(n, a, h); 
    Vect p = solver::GetLineNearest(x, n, a, h);
    std::cerr 
        << "n=" << n
        << " u=" << u
        << " h=" << h
        << " a=" << a
        << " x0=" << xx[0]
        << " x1=" << xx[1]
        << " p=" << p
        << " pe=" << pe
        << std::endl;
    assert((p - pe).norm() < 1e-8);
  };

  Scal nx, ny, u, hx, hy;
  Scal nx0, ny0, u0, hx0, hy0;

  std::cerr << "check nearest point to line" << std::endl;
  nx0 = 1.; ny0 = 0.; u0 = 0.5; hx0 = 1.; hy0 = 1.; 
  nx = nx0; ny = ny0; u = u0; hx = hx0; hy = hy0;

  std::cerr << "base" << std::endl;
  f(Vect(0., 0., 0.), nx, ny, u, hx, hy, Vect(0., 0., 0.));
  f(Vect(1., 0., 0.), nx, ny, u, hx, hy, Vect(0., 0., 0.));
  f(Vect(1., 0., 0.), nx, ny, 0.75, hx, hy, Vect(0.25, 0., 0.));
  f(Vect(1., 1., 0.), nx, ny, u, hx, hy, Vect(0., 0.5, 0.));
  f(Vect(1.1, 0.9, 0.), 1., 1., u, hx, hy, Vect(0.1, -0.1, 0.));
}

void Plot() {
  using namespace solver;
  Vect n(0.1, 0.2, 0.29);
  {
    std::string fn = "vof_a.dat";
    std::ofstream f(fn);
    for (Scal u = 0.; u < 1.+1e-6; u += 0.001) {
      f << u << " " << GetLineA1(n, u) << "\n";
    }
  }

  {
    std::string fn = "vof_u.dat";
    std::ofstream f(fn);
    Scal nh = n.norm1() * 0.5;
    for (Scal a = -nh; a < nh+1e-6; a += 0.001 * nh * 2.) {
      f << a << " " << GetLineU1(n, a) << "\n";
    }
  }
}

void TestFit() {
  std::cerr << "check plane fitting" << std::endl;
  using namespace solver;

  using V = Vect;
  using VV = GVect<V, 3>;

  auto f = [](const VV& a, const V& xe) {
    auto x = SolveSingular(a);
    std::cout 
        << "a=" << a
        << " x=" << x 
        << " xe=" << xe
        << std::endl;
    assert(x.dist(xe) < 1e-12);
  };

  f(VV(V(-2., 1., 1.), V(1., -2., 1.), V(1., 1., -2.)), V(1., 1., 1.));
  f(VV(V(1., 0., 0.), V(0., 1., 0.), V(1., 1., 0.)), V(0., 0., 1.));

  auto g = [](const std::vector<V>& xx, const V& ne) {
    auto n = GetFitN(xx);
    std::cout 
        << "xx=" << xx
        << " n=" << n 
        << " ne=" << ne
        << std::endl;
    assert(n.dist(ne) < 1e-12);
  };

  g({V(0., 0., 0.), V(0., 1., 0.), V(0., 1., 1.), V(0., 0., 1.)}, V(1., 0., 0.));
}

void TestNearest() {
  std::cerr << "check nearest to polygon" << std::endl;
  using namespace solver;
  using V = Vect;

  auto f1 = [](V x, V x0, V x1, V n, V xne) {
    auto xn = GetNearestHalf(x, x0, x1, n);
    std::cout << "x0x1:"
        << " x=" << x 
        << " x0=" << x0
        << " x1=" << x1
        << " n=" << n
        << " xn=" << xn
        << " xne=" << xne
        << std::endl;
    assert(xne.dist(xn) < 1e-12);
  };

  auto f2 = [](V x, V x0, V x1, V xh, V n, V xne) {
    auto xn = GetNearestHalf(x, x0, x1, xh, n);
    std::cout << "xh:"
        << " x=" << x 
        << " x0=" << x0
        << " x1=" << x1
        << " xh=" << xh
        << " n=" << n
        << " xn=" << xn
        << " xne=" << xne
        << std::endl;
    assert(xne.dist(xn) < 1e-12);
  };

  auto f3 = [](V x, const std::vector<V>& xx, V n, V xne) {
    auto xn = GetNearest(x, xx, n);
    std::cout << "tri:"
        << " x=" << x 
        << " xx=" << xx
        << " n=" << n
        << " xn=" << xn
        << " xne=" << xne
        << std::endl;
    assert(xne.dist(xn) < 1e-12);
  };

  f1(V(1.,1.,7.), V(-8.,0.,0.), V(7.,0.,0.), V(0.,0.,3.), V(1.,1.,0.));
  f1(V(1.,1.,9.), V(12.,0.,0.), V(-3.,0.,0.), V(0.,0.,7.), V(1.,0.,0.));
  f1(V(1.,1.,7.), V(0.,0.,0.), V(0.5,0.,0.), V(0.,0.,3.), V(0.5,1.,0.));
  f1(V(-1.,1.,7.), V(0.,0.,0.), V(0.5,0.,0.), V(0.,0.,3.), V(0.,1.,0.));
  f1(V(1.,1.,7.), V(0.5,0.,0.), V(0.,0.,0.), V(0.,0.,3.), V(0.5,0.,0.));

  f2(V(1.,1.,7.), V(-1.,0.,0.), V(11.,0.,0.), V(0.,1.,0.), V(0.,0.,3.), V(1.,1.,0.));
  f2(V(1.,1.,7.), V(-2.,0.,0.), V(11.,0.,0.), V(0.,-1.,0.), V(0.,0.,3.), V(1.,0.,0.));

  {
    std::vector<V> xx = {V(0.,0.,0.), V(1.,0.,0.), V(0.,1.,0.)};
    V n(0.,0.,1.);
    f3(V(10.,10.,1.), xx, n, V(0.5,0.5,0.));
    f3(V(10.,0.,1.), xx, n, V(1.,0.,0.));
    f3(V(0.1,0.2,1.), xx, n, V(0.1,0.2,0.));
  }
}


int main() {
  TestNearest();
  return 0;
  TestFit();
  Plot();
  TestRect();
  TestRandom();
  TestVol();
  TestVolStr();
  TestCenter();
  TestLineNearest();
}
