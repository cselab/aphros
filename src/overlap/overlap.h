// Created by Petr Karnakov on 02.04.2019
// Copyright 2019 ETH Zurich

#include "geom/vect.h"

// Returns volume of intersection of cell and sphere.
// x: cell center
// h: cell size
// c: sphere center
// r: radius
double GetSphereOverlap(
    const generic::Vect<double, 3>& x, const generic::Vect<double, 3>& h,
    const generic::Vect<double, 3>& c, double r);
