// Created by Sergey Litvinov on 24.12.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <aphros/march/march.h>

int main() {
  enum { X, Y, Z };
  enum { N = 3 * 3 * MARCH_NTRI };
  double u[] = {-1, 1, 1, 1, 1, 1, 1, -1}, tri[N], offset[N], x, a, b;
  int first[N], second[N], i, n;

  i = 0;
  march_cube_location(u, &n, tri, first, second, offset);
  a = MARCH_O[first[i]][X];
  b = MARCH_O[second[i]][X];
  x = a + (b - a) * offset[i];
  assert(fabs(x - tri[3 * i + X]) < 1e-12);
}
