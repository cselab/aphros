#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include "mc.h"

#define T2                                                    \
  (((x * x + y * y + z * z + R * R - r * r) *                 \
        (x * x + y * y + z * z + R * R - r * r) -             \
    4 * R * R * (x * x + y * y)) *                            \
   ((x * x + (y + R) * (y + R) + z * z + R * R - r * r) *     \
        (x * x + (y + R) * (y + R) + z * z + R * R - r * r) - \
    4 * R * R * ((y + R) * (y + R) + z * z)))
#define MC                                                                    \
  (-26.5298 * (1 - x) * (1 - y) * (1 - z) + 81.9199 * x * (1 - y) * (1 - z) - \
   100.68 * x * y * (1 - z) + 3.5498 * (1 - x) * y * (1 - z) +                \
   24.1201 * (1 - x) * (1 - y) * z - 74.4702 * x * (1 - y) * z +              \
   91.5298 * x * y * z - 3.22998 * (1 - x) * y * z)
#define CHAIR                                    \
  (x * x + y * y + z * z - 0.95f * 25) *         \
          (x * x + y * y + z * z - 0.95f * 25) - \
      0.8f * ((z - 5) * (z - 5) - 2 * x * x) * ((z + 5) * (z + 5) - 2 * y * y)

#define CUSHIN                                                             \
  (z * z * x * x - z * z * z * z - 2 * z * x * x + 2 * z * z * z + x * x - \
   z * z - (x * x - z) * (x * x - z) - y * y * y * y - 2 * x * x * y * y - \
   y * y * z * z + 2 * y * y * z + y * y)

enum { X, Y, Z };
static void obj(int, double*, int, int*);

int main() {
  int nx = 40, ny = 40, nz = 40;
  int i, j, k, m;
  double x, y, z;
  double sx, sy, sz;
  double tx, ty, tz;
  double r, R;
  double* buf;
  int nv, nt;
  double* ver;
  int* tri;
  struct MarchLewiner* march;

  buf = (double*)malloc(nx * ny * nz * sizeof(*buf));
  if (buf == NULL) {
    fprintf(stderr, "alloc failed\n");
    exit(2);
  }
  r = 1.85;
  R = 4;
  sx = nx / 16.0;
  sy = ny / 16.0;
  sz = nz / 16.0;
  tx = nx / (2 * sx);
  ty = ny / (2 * sy) + 1.5;
  tz = nz / (2 * sz);
  m = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++) {
        x = i / sx - tx;
        y = j / sy - ty;
        z = k / sz - tz;
        buf[m++] = T2;
      }

  march = march_lewiner_ini(nx, ny, nz);
  march_lewiner_apply(march, buf, &nv, &ver, &nt, &tri);
  obj(nv, ver, nt, tri);
  march_lewiner_fin(march);
  free(buf);
  return 0;
}

static void obj(int nv, double* ver, int nt, int* tri) {
  int i;

  printf("# File type: ASCII OBJ\n");
  for (i = 0; i < nv; i++)
    printf(
        "v %.16g %.16g %.16g\n", ver[3 * i + X], ver[3 * i + Y],
        ver[3 * i + Z]);
  for (i = 0; i < nt; i++)
    printf(
        "f %d %d %d\n", tri[3 * i + X] + 1, tri[3 * i + Y] + 1,
        tri[3 * i + Z] + 1);
}
