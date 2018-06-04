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

void TestRandom() {
  std::cerr << "check u=u(a(u)) on random input" << std::endl;
  std::default_random_engine g(0);
  std::uniform_real_distribution<double> f(0.,1.);
  size_t ni = 10000;
  Scal e = 0.; // error
  for (size_t i = 0; i < ni; ++i) {
    while (true) {
      Vect n(f(g), f(g), 0.);
      n /= n.norm();
      Scal u = f(g);
      e += std::abs(solver::GetLineU1(n, solver::GetLineA1(n, u)) - u);

      Vect h(0.1, 0.2, 0.3);
      e += std::abs(solver::GetLineU(n, solver::GetLineA(n, u, h), h) - u);
      break;
    }
  }
  e /= ni;
  std::cerr << "TestLine(): error=" << e << " samples=" << ni << std::endl;
  assert(e < 1e-16);
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
    n /= n.norm();
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
  v(1., 1., 1., -Vect(0.5).norm(), 0.);
  v(1., 1., 1., Vect(0.5).norm(), 1.);
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
    if (xx) {
      dv = solver::GetLineVolX(n, a, h, dx);
    } else {
      dv = solver::GetLineVolY(n, a, h, dx);
    }
    std::cerr 
        << "n=" << n
        << " u=" << u
        << " h=" << h
        << " dx=" << dx
        << " a=" << a
        << " dv=" << dv
        << " dve=" << dve
        << std::endl;
    assert(std::abs(dv - dve) < 1e-12);
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
    if (xx) {
      dv = solver::GetLineVolStrX(n, a, h, dx, dxu);
    } else {
      dv = solver::GetLineVolStrY(n, a, h, dx, dxu);
    }
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
    assert(std::abs(dv - dve) < 1e-12);
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
    assert((c - ce).norm() < 1e-12);
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
void TestNearest() {
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
    Vect p = solver::GetNearest(x, n, a, h);
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
    assert((p - pe).norm() < 1e-12);
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

int main() {
  TestRect();
  return 0;
  TestRect();
  TestRandom();
  TestVol();
  TestVolStr();
  TestCenter();
  TestNearest();
}
