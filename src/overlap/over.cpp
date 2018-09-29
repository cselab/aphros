#include "overlap.hpp"

#include "over.h"

double GetSphereOverlap(const GVect<double, 3>& x, const GVect<double, 3>& h, 
                        const GVect<double, 3>& c, double r) {
  std::array<vector_t, 8> vv{
    vector_t{-1, -1, -1},
    vector_t{ 1, -1, -1},
    vector_t{ 1,  1, -1},
    vector_t{-1,  1, -1},
    vector_t{-1, -1,  1},
    vector_t{ 1, -1,  1},
    vector_t{ 1,  1,  1},
    vector_t{-1,  1,  1}
  };

  auto hh = h * 0.5;
  for (size_t i = 0; i < vv.size(); ++i) {
    auto& v = vv[i];
    v[0] = x[0] + v[0] * hh[0];
    v[1] = x[1] + v[1] * hh[1];
    v[2] = x[2] + v[2] * hh[2];
  }

  Hexahedron hex{vv};
  Sphere s{vector_t{c[0], c[1], c[2]}, r};

  return overlap(s, hex) / h.prod();
}
