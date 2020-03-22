// Created by Petr Karnakov on 04.03.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include <util/posthook.h>

template <class M>
void InitVelHook(
    FieldCell<typename M::Vect>& fcvel, const Vars&, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  if (m.IsRoot()) {
    std::cout << "initial velocity: Bell,Colella,Glaz 1989" << std::endl;
  }
  for (auto c : m.AllCells()) {
    // psi = 1/pi * sin**2(pi*x) * sin**2(pi*y)
    const Vect xx = m.GetCenter(c);
    const Scal x = xx[0];
    const Scal y = xx[1];
    Scal& u = fcvel[c][0];
    Scal& v = fcvel[c][1];
    const Scal pi = M_PI;
    using std::sin;
    using std::cos;
    u = -2 * sqr(sin(pi * x)) * sin(pi * y) * cos(pi * y);
    v = 2 * sqr(sin(pi * y)) * sin(pi * x) * cos(pi * x);
  }
}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

template void InitVelHook(FieldCell<Vect>&, const Vars&, const M&);
