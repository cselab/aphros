#include <iostream>

#include "overlap.hpp"

int main() {
  while (true) {
    scalar_t cx, cy, cz, h;

    std::cin >> cx >> cy >> cz >> h;

    if (std::cin.fail()) {
      break;
    }

    
    scalar_t a = h * 0.5;

    scalar_t x0, x1, y0, y1, z0, z1;

    x0 = cx - a; x1 = cx + a;
    y0 = cy - a; y1 = cy + a;
    z0 = cz - a; z1 = cz + a;

    vector_t v0{x0, y0, z0};
    vector_t v1{x1, y0, z0};
    vector_t v2{x1, y1, z0};
    vector_t v3{x0, y1, z0};
    vector_t v4{x0, y0, z1};
    vector_t v5{x1, y0, z1};
    vector_t v6{x1, y1, z1};
    vector_t v7{x0, y1, z1};

    Hexahedron hex{v0, v1, v2, v3, v4, v5, v6, v7};

    Sphere s{vector_t::Constant(0), 1};

    scalar_t r = overlap(s, hex);

    std::cout.precision(20);
    std::cout << r << std::endl;
  }
}
