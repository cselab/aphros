#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <march.h>

enum { X, Y, Z };
enum {N = 3 * 3 * MARCH_NTRI};
static double av(double a, double b, double o) {
  return a + (b - a) * o;
}

int main() {
  double u[] = {-1, 1, 1, 1, 1, 1, 1, -1}, tri[N], offset[N], x;
  int first[N], second[N], i, n;

  i = 0;
  march_cube_location(u, &n, tri, first, second, offset);
  x = av(MARCH_O[first[i]][X], MARCH_O[second[i]][X], offset[i]);
  assert(fabs(x - tri[3 * i + X] < 1e-12));
}
