#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <march.h>

enum { X, Y, Z };
enum {N = MARCH_NTRI};
static double av(double a, double b, double o) {
  return a + (b - a) * o;
}

int main() {
  double u[] = {-1, 0, 0, 0, 0, 0, 0, 1}, tri[3 * 3 * N], offset[3 * N], x;
  int first[3 * N], second[3 * N], i, n;

  i = 0;
  march_cube_location(u, &n, tri, first, second, offset);
  x = av(MARCH_O[first[i]][X], MARCH_O[second[i]][X], offset[i]);
  assert(fabs(x - tri[3 * i + X] < 1e-12));
}
