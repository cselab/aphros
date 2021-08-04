#include <float.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <array>
#include <bitset>
#include <cassert>
#include <numeric>
#include <type_traits>

#include "double_prec.inc"
#include "util.inc"
#include "overlap.inc"

int main() {
  vector_t a{1, 2, 3};
  vector_t b{5, 8, 3};
  vector_t c;

  c = gramSchmidt(a, b);
  printf("%.16e %.16e %.16e\n", c[0], c[1], c[2]);
}
