#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

enum { N = 1024 };

#define MALLOC(n, p)                                                           \
  do {                                                                         \
    *(p) = malloc(n * sizeof(**(p)));                                          \
    if (*(p) == NULL) {                                                        \
      fprintf(stderr, "%s:%d: alloc failed, n = %d\n", __FILE__, __LINE__, n); \
      exit(2);                                                                 \
    }                                                                          \
  } while (0)

struct Mesh {
  int nt, nv;
  float* r;
  int* t;
};

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

static int tri_center(
    float* a, float* b, float* c, /**/ double* x, double* y, double* z) {
  enum { X, Y, Z };
  *x += (a[X] + b[X] + c[X]) / 3;
  *y += (a[Y] + b[Y] + c[Y]) / 3;
  *z += (a[Z] + b[Z] + c[Z]) / 3;
  return 0;
}

static double tri_dot(
    float* a, float* b, float* c, double x, double y, double z) {
  enum { X, Y, Z };
  double bx, by, bz, cx, cy, cz, nx, ny, nz;

  x -= a[X];
  y -= a[Y];
  z -= a[Z];
  bx = b[X] - a[X];
  by = b[Y] - a[Y];
  bz = b[Z] - a[Z];
  cx = c[X] - a[X];
  cy = c[Y] - a[Y];
  cz = c[Z] - a[Z];

  nx = by * cz - bz * cy;
  ny = bz * cx - bx * cz;
  nz = bx * cy - by * cx;

  x -= (bx + cx) / 3;
  y -= (by + cy) / 3;
  z -= (bz + cz) / 3;

  return nx * x + ny * y + nz * z;
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
  if ((int)fread(c, sizeof(*c), n, f) != n) {
    fprintf(stderr, "%s:%d: failt to read, n = %d\n", __FILE__, __LINE__, n);
    exit(2);
  }
  return 0;
}

static int vector(FILE* f, int n, char* name, float* c) {
  char s[N];
  line(s, f);
  if (sscanf(s, "VECTORS %s float", name) != 1) {
    fprintf(stderr, "%s:%d: expect VECTOR, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  if ((int)fread(c, sizeof(*c), 3 * n, f) != 3 * n) {
    fprintf(stderr, "%s:%d: failt to read, n = %d\n", __FILE__, __LINE__, n);
    exit(2);
  }
  return 0;
}

static int get1(float* r, int i, /**/ float* a) {
  enum { X, Y, Z };
  a[X] = r[3 * i + X];
  a[Y] = r[3 * i + Y];
  a[Z] = r[3 * i + Z];
  return 0;
}

static int get3(struct Mesh* mesh, int t, /**/ float* a, float* b, float* c) {
  int nt, *tri, i, j, k;
  float* r;
  nt = mesh->nt;
  tri = mesh->t;
  r = mesh->r;
  if (t > nt) {
    fprintf(stderr, "%s:%d: t=%d > nt=%d\n", __FILE__, __LINE__, t, nt);
    exit(1);
  }

  i = tri[3 * t];
  j = tri[3 * t + 1];
  k = tri[3 * t + 2];
  get1(r, i, a);
  get1(r, j, b);
  get1(r, k, c);
  return 0;
}

static int* u_root;
static int max_arg(int n, int* a) {
  int i, j, m;
  j = 0;
  m = a[j];
  for (i = 0; i < n; i++) {
    if (a[i] > m) {
      m = a[i];
      j = i;
    }
  }
  return j;
}

static void u_ini(int n) {
  MALLOC(n, &u_root);
  while (n--)
    u_root[n] = n;
}

static void u_fin() {
  free(u_root);
}

static int u_find(int v) {
  if (v == u_root[v]) return v;
  return u_root[v] = u_find(u_root[v]);
}

static void u_union(int a, int b) {
  a = u_find(a);
  b = u_find(b);
  if (a != b) u_root[b] = a;
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

int main() {
  FILE* f;
  char s[N], name[N];
  int nv, nt, nb, u, v, w;
  float *r, *nn;
  double *cx, *cy, *cz, *dot;
  int *t, *t0, *a, *b, *id, *c, *cnt;
  int i, j, k;
  float x[3], y[3], z[3], *cl;
  struct Mesh mesh;

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
  if ((int)fread(r, sizeof(*r), 3 * nv, f) != 3 * nv) {
    fprintf(stderr, "%s:%d: failt to read, nv = %d\n", __FILE__, __LINE__, nv);
    exit(2);
  }
  swap(3 * nv, sizeof(*r), r);

  while (line(s, f) == 0 && s[0] == '\0')
    ;
  if (sscanf(s, "POLYGONS %d %*d", &nt) != 1) {
    fprintf(
        stderr, "%s:%d: expect POLYGONS, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  MALLOC(4 * nt, &t0);
  MALLOC(3 * nt, &t);
  if ((int)fread(t0, sizeof(*t0), 4 * nt, f) != 4 * nt) {
    fprintf(stderr, "%s:%d: failt to read, nt = %d\n", __FILE__, __LINE__, nt);
    exit(2);
  }
  for (i = 0, a = t, b = t0; i < nt; i++) {
    b++;
    *a++ = *b++;
    *a++ = *b++;
    *a++ = *b++;
  }
  swap(3 * nt, sizeof(*t), t);
  while (line(s, f) == 0 && s[0] == '\0')
    ;
  if (sscanf(s, "CELL_DATA %*d") != 0) {
    fprintf(
        stderr, "%s:%d: expect CELL_DATA, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  MALLOC(nt, &cl);
  MALLOC(nt, &c);
  for (;;) {
    scalar(f, nt, name, cl);
    if (eq(name, "cl")) break;
  }
  swap(nt, sizeof(*cl), cl);

  while (line(s, f) == 0 && s[0] == '\0')
    ;
  if (sscanf(s, "POINT_DATA %*d") != 0) {
    fprintf(
        stderr, "%s:%d: expect POINT_DATA, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  MALLOC(3 * nv, &nn);
  for (;;) {
    vector(f, nv, name, nn);
    if (eq(name, "nn")) break;
  }

  u_ini(nv);
  for (i = 0; i < nt; i++) {
    if (cl[i] < 0) continue;
    u = t[3 * i];
    v = t[3 * i + 1];
    w = t[3 * i + 2];
    u_union(u, v);
    u_union(v, w);
    u_union(u, w);
  }
  MALLOC(nt, &id);
  for (i = 0; i < nt; i++) {
    u = t[3 * i];
    c[i] = u_find(u);
  }
  for (i = 0; i < nt; i++)
    id[i] = -1;
  for (j = i = 0; i < nt; i++) {
    k = c[i];
    if (id[k] == -1) {
      id[k] = j++;
    }
    c[i] = id[k];
  }
  nb = j;
  MALLOC(nb, &cnt);
  for (i = 0; i < nb; i++)
    cnt[i] = 0;
  for (i = 0; i < nt; i++)
    cnt[c[i]]++;

  k = max_arg(nb, cnt); /* water = 0 */
  for (i = 0; i < nt; i++) {
    if (c[i] == k)
      c[i] = 0;
    else if (c[i] == 0)
      c[i] = k;
  }

  mesh.r = r;
  mesh.t = t;
  mesh.nt = nt;
  mesh.nv = nv;
  MALLOC(nb, &cx);
  MALLOC(nb, &cy);
  MALLOC(nb, &cz);
  MALLOC(nb, &dot);
  for (i = 0; i < nb; i++)
    cx[i] = cy[i] = cz[i] = dot[i] = 0;
  for (i = 0; i < nt; i++) {
    k = c[i];
    get3(&mesh, i, x, y, z);
    tri_center(x, y, z, &cx[k], &cy[k], &cz[k]);
  }
  for (i = 0; i < nb; i++) {
    cx[i] /= cnt[i];
    cy[i] /= cnt[i];
    cz[i] /= cnt[i];
  }
  for (i = 0; i < nt; i++) {
    k = c[i];
    get3(&mesh, i, x, y, z);
    dot[k] += tri_dot(x, y, z, cx[k], cy[k], cz[k]);
  }

  for (i = 0; i < nt; i++) {
    k = c[i];
    if (k != 0 && dot[k] > 0) c[i] = 0;
  }

  for (i = 0; i < nt; i++)
    if (cl[i] < 0) c[i] = -1;

  printf(
      "# vtk DataFile Version 2.0\n"
      "Interface from marching cubes\n"
      "BINARY\n"
      "DATASET POLYDATA\n"
      "POINTS %d float\n",
      nv);
  swap(3 * nv, sizeof(*r), r);
  fwrite(r, sizeof(*r), 3 * nv, stdout);
  printf("POLYGONS %d %d\n", nt, 4 * nt);
  fwrite(t0, sizeof(*t0), 4 * nt, stdout);
  printf(
      "CELL_DATA %d\n"
      "SCALARS cl int\n"
      "LOOKUP_TABLE default\n",
      nt);
  swap(nt, sizeof(*c), c);
  fwrite(c, sizeof(*c), nt, stdout);
  printf(
      "POINT_DATA %d\n"
      "VECTORS nn float\n",
      nv);
  fwrite(nn, sizeof(*nn), 3 * nv, stdout);

  u_fin();
  free(c);
  free(cl);
  free(cnt);
  free(id);
  free(cx);
  free(cy);
  free(cz);
  free(dot);
  free(nn);
  free(r);
  free(t);
  free(t0);
}
