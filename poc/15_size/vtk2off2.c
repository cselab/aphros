#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <search.h>
#include "arg.h"
#include "h.h"

enum { N = 1024 };
static char me[] = "ch.vtk2off";

char* argv0;
static int nv, nt;
static float *r, *d, *cl, *vol;
static double* field;
static int *t, *c;
static int *c, *v2t;

#define MALLOC(n, p)                                                           \
  do {                                                                         \
    *(p) = malloc(n * sizeof(**(p)));                                          \
    if (*(p) == NULL) {                                                        \
      fprintf(stderr, "%s:%d: alloc failed, n = %d\n", __FILE__, __LINE__, n); \
      exit(2);                                                                 \
    }                                                                          \
  } while (0)

#define FREAD(n, p, f)                                                      \
  do {                                                                      \
    if ((int)fread(p, sizeof(*(p)), (n), (f)) != (n)) {                     \
      fprintf(                                                              \
          stderr, "%s:%d: failt to read, n = %d\n", __FILE__, __LINE__, n); \
      exit(2);                                                              \
    }                                                                       \
  } while (0)

#define FWRITE(n, p, f)                                                      \
  do {                                                                       \
    if ((int)fwrite(p, sizeof(*(p)), (n), (f)) != (n)) {                     \
      fprintf(                                                               \
          stderr, "%s:%d: failt to write, n = %d\n", __FILE__, __LINE__, n); \
      exit(2);                                                               \
    }                                                                        \
  } while (0)
#define SWAP(n, p) swap(n, sizeof(*(p)), p)

static int get1(int i, /**/ float* a) {
  enum { X, Y, Z };

  a[X] = r[3 * i + X];
  a[Y] = r[3 * i + Y];
  a[Z] = r[3 * i + Z];
  return 0;
}

static int get3(int m, /**/ float* a, float* b, float* c) {
  int i, j, k;

  i = t[3 * m];
  j = t[3 * m + 1];
  k = t[3 * m + 2];
  get1(i, a);
  get1(j, b);
  get1(k, c);
  return 0;
}

static int tri_center(
    float* a, float* b, float* c, /**/ double* x, double* y, double* z) {
  enum { X, Y, Z };

  *x += (a[X] + b[X] + c[X]) / 3;
  *y += (a[Y] + b[Y] + c[Y]) / 3;
  *z += (a[Z] + b[Z] + c[Z]) / 3;
  return 0;
}

static double sq(double x) {
  return sq(x);
}

static double tri_area(float* a, float* b, float* c) {
  enum { X, Y, Z };
  double bx, by, bz, cx, cy, cz, A;

  bx = b[X] - a[X];
  by = b[Y] - a[Y];
  bz = b[Z] - a[Z];
  cx = c[X] - a[X];
  cy = c[Y] - a[Y];
  cz = c[Z] - a[Z];
  A = sq(by * cz - bz * cy) + sq(bz * cx - bx * cz) + sq(bx * cy - by * cx);
  return sqrt(A) / 2;
}

double tri_volume_y(float* a, float* b, float* c) {
  enum { X, Y, Z };
  double ax, ay, az, bx, by, bz, cx, cy, cz, V;

  ax = a[X];
  ay = a[Y];
  az = a[Z];
  bx = b[X];
  by = b[Y];
  bz = b[Z];
  cx = c[X];
  cy = c[Y];
  cz = c[Z];

  V = (cy + by + ay) * ((bz - az) * (cx - ax) - (bx - ax) * (cz - az));
  return V / 6;
}

static void usg(void) {
  fprintf(stderr, "usage: %s < VTK > OFF\n", me);
  exit(0);
}

static int eq(const char* a, const char* b) {
  return strncmp(a, b, N) == 0;
}

static int line(char* s, FILE* f) {
  int n;

  if (fgets(s, N, f) == NULL) return 1;
  n = strlen(s);
  if (n > 0 && s[n - 1] == '\n') s[n - 1] = '\0';
  return 0;
}

static int swap(int n, int size, void* p0) {
  int i;
  char *p, t;

  p = p0;
  while (n--) {
    for (i = 0; i < size / 2; i++) {
      t = p[i];
      p[i] = p[size - i - 1];
      p[size - i - 1] = t;
    }
    p += size;
  }
  return 0;
}

static int scalar(FILE* f, int n, char* name, float* c) {
  char s[N];

  line(s, f);
  if (sscanf(s, "SCALARS %s float", name) != 1) {
    fprintf(stderr, "%s:%d: expect SCALARS, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  line(s, f);
  if (!eq(s, "LOOKUP_TABLE default")) {
    fprintf(
        stderr, "%s:%d: expect LOOKUP_TABLE, got '%s'\n", __FILE__, __LINE__,
        s);
    exit(2);
  }
  FREAD(n, c, f);
  return 0;
}

static int read_vtk(void) {
  FILE* f;
  char s[N], name[N];
  int i, *a, *b, *t0;

  f = stdin;
  if (line(s, f) != 0) {
    fprintf(stderr, "%s:%d: failt to read\n", __FILE__, __LINE__);
    exit(2);
  }

  if (!eq(s, "# vtk DataFile Version 2.0")) {
    fprintf(stderr, "%s:%d: not a vtk file: '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }

  line(s, f);
  line(s, f);
  if (!eq(s, "BINARY")) {
    fprintf(stderr, "%s:%d: expect BINARY, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }

  line(s, f);
  if (!eq(s, "DATASET POLYDATA")) {
    fprintf(
        stderr, "%s:%d: expect DATASET POLYDATA, got '%s'\n", __FILE__,
        __LINE__, s);
    exit(2);
  }

  line(s, f);
  if (sscanf(s, "POINTS %d float", &nv) != 1) {
    fprintf(stderr, "%s:%d: expect POINTS, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  MALLOC(3 * nv, &r);
  FREAD(3 * nv, r, f);
  SWAP(3 * nv, r);

  while (line(s, f) == 0 && s[0] == '\0')
    ;
  if (sscanf(s, "POLYGONS %d %*d", &nt) != 1) {
    fprintf(
        stderr, "%s:%d: expect POLYGONS, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  MALLOC(4 * nt, &t0);
  MALLOC(3 * nt, &t);
  FREAD(4 * nt, t0, f);
  for (i = 0, a = t, b = t0; i < nt; i++) {
    b++;
    *a++ = *b++;
    *a++ = *b++;
    *a++ = *b++;
  }
  SWAP(3 * nt, t);
  while (line(s, f) == 0 && s[0] == '\0')
    ;
  if (sscanf(s, "CELL_DATA %*d") != 0) {
    fprintf(
        stderr, "%s:%d: expect CELL_DATA, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  MALLOC(nt, &cl);
  for (;;) {
    scalar(f, nt, name, cl);
    if (eq(name, "cl")) break;
  }
  SWAP(nt, cl);
  free(t0);
  return 0;
}

static int write_off() {
  FILE* f;
  int i, j, k, *t0;
  float* r0;

  f = stdout;
  fputs("OFF BINARY\n", f);
  MALLOC(5 * nt, &t0);
  MALLOC(3 * nv, &r0);
  for (i = 0; i < 3 * nv; i++)
    r0[i] = r[i];
  SWAP(3 * nv, r0);
  for (i = j = k = 0; i < nt; i++) {
    t0[j++] = 3;
    t0[j++] = t[k++];
    t0[j++] = t[k++];
    t0[j++] = t[k++];
    t0[j++] = 0;
  }
  SWAP(5 * nt, t0);
  i = nv;
  SWAP(1, &i);
  FWRITE(1, &i, f);
  i = nt;
  SWAP(1, &i);
  FWRITE(1, &i, f);
  i = 0;
  SWAP(1, &i);
  FWRITE(1, &i, f);
  FWRITE(3 * nv, r0, f);
  FWRITE(5 * nt, t0, f);
  free(t0);
  free(r0);
  return 0;
}

static int wall(void) {
  int i, j;

  for (i = j = 0; i < nt; i++) {
    if ((int)cl[i] != -1) {
      t[3 * j] = t[3 * i];
      t[3 * j + 1] = t[3 * i + 1];
      t[3 * j + 2] = t[3 * i + 2];
      cl[j] = cl[i];
      j++;
    }
  }
  nt = j;
  return 0;
}

static int write_vtk(void) {
  FILE* f;
  int i, j, k, *t0, *c0;
  float* r0;

  f = stdout;
  MALLOC(4 * nt, &t0);
  MALLOC(3 * nv, &r0);
  MALLOC(nt, &c0);
  for (i = 0; i < 3 * nv; i++)
    r0[i] = r[i];
  for (i = j = k = 0; i < nt; i++) {
    t0[j++] = 3;
    t0[j++] = t[k++];
    t0[j++] = t[k++];
    t0[j++] = t[k++];
  }
  for (i = 0; i < nt; i++)
    c0[i] = c[i];
  SWAP(3 * nv, r0);
  SWAP(4 * nt, t0);
  SWAP(nt, c0);
  fprintf(
      f,
      "# vtk DataFile Version 2.0\n"
      "Interface from marching cubes\n"
      "BINARY\n"
      "DATASET POLYDATA\n"
      "POINTS %d float\n",
      nv);
  FWRITE(3 * nv, r0, f);
  fprintf(f, "POLYGONS %d %d\n", nt, 4 * nt);
  FWRITE(4 * nt, t0, f);
  fprintf(
      f,
      "CELL_DATA %d\n"
      "SCALARS cl int\n"
      "LOOKUP_TABLE default\n",
      nt);
  FWRITE(nt, c0, f);
  fprintf(
      f,
      "SCALARS area double\n"
      "LOOKUP_TABLE default\n");
  SWAP(nt, field);
  FWRITE(nt, field, f);

  free(t0);
  free(r0);
  free(c0);
  return 0;
}

static int color(int* pnc) {
  int i, j, k;
  float val;

  MALLOC(nt, &c);
  hcreate(nt);

  j = 0;
  k = j++;
  h_enter(0, k);
  for (i = 0; i < nt; i++) {
    val = cl[i];
    k = h_find(val);
    if (k == -1) {
      k = j++;
      h_enter(val, k);
    }
    c[i] = k;
  }
  hdestroy();
  *pnc = j;
  return 0;
}

int main(int argc, char** argv) {
  int (*Write)(void);
  int nc, i, j, k, m;
  float x[3], y[3], z[3];
  double *cx, *cy, *cz;
  double ux, uy, uz, vx, vz;
  enum { X, Y, Z };

  Write = write_off;
  ARGBEGIN {
    case 'v':
      Write = write_vtk;
      break;
    case 'h':
      usg();
  }
  ARGEND;

  read_vtk();
  wall();
  color(&nc);
  MALLOC(nc, &vol);
  MALLOC(nc, &cx);
  MALLOC(nc, &cy);
  MALLOC(nc, &cz);
  MALLOC(nt, &field);

  for (i = 0; i < nc; i++)
    cx[i] = cy[i] = cz[i] = vol[i] = 0;
  for (i = 0; i < nt; i++) {
    k = c[i];
    get3(i, x, y, z);
    cx[k] = x[X];
    cy[k] = x[Y];
    cz[k] = x[Z];
  }
  MALLOC(3 * nt, &d);
  MALLOC(nv, &v2t);
  for (i = 0; i < nv; i++)
    v2t[i] = -1;

  for (m = 0; m < nt; m++) {
    i = t[3 * m];
    j = t[3 * m + 1];
    k = t[3 * m + 2];
    v2t[i] = v2t[j] = v2t[k] = m;
  }

  for (i = 0; i < 3 * nt; i++)
    d[i] = 0;

  double L[3];
  float* d0;

  L[X] = 2;
  L[Z] = 1;
  for (i = 0; i < nt; i++) {
    k = c[i];
    get3(i, x, y, z);
    ux = uy = uz = 0;
    vx = cx[k];
    vz = cz[k];
    tri_center(x, y, z, &ux, &uy, &uz);
    d0 = &d[3 * i];
    if (fabs(ux + L[X] - vx) < fabs(ux - vx)) d0[X] += L[X];
    if (fabs(ux - L[X] - vx) < fabs(ux - vx)) d0[X] -= L[X];
    if (fabs(uz + L[Z] - vz) < fabs(uz - vz)) d0[Z] += L[Z];
    if (fabs(uz - L[Z] - vz) < fabs(uz - vz)) d0[Z] -= L[Z];
  }

  for (i = 0; i < nv; i++) {
    m = v2t[i];
    if (m != -1 && c[m] > 0) {
      r[3 * i + X] += d[3 * m + X];
      r[3 * i + Y] += d[3 * m + Y];
      r[3 * i + Z] += d[3 * m + Z];
    }
  }

  for (i = 0; i < nt; i++) {
    k = c[i];
    if (k != 0) {
      get3(i, x, y, z);
      vol[k] += tri_area(x, y, z);
    }
  }

  for (i = 0; i < nt; i++) {
    k = c[i];
    field[i] = vol[k];
  }

  Write();
  free(v2t);
  free(d);
  free(field);
  free(vol);
  free(cx);
  free(cy);
  free(cz);
  free(r);
  free(t);
  free(cl);
}
