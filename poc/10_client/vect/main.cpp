// Created by Sergey Litvinov on 11.03.2019
// Copyright 2019 ETH Zurich

#include "geom/mesh.h"

int main() {
  using std::cout;
  GVect<double, 3> a(1), b(2), c;

  c = a + b;

  cout << a << '\n';
  cout << b << '\n';
  cout << c << '\n';

  cout << a.dot(b) << '\n';
  cout << a.cross(b) << '\n';

  return 0;
}
