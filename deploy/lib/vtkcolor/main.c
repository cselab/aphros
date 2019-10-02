#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {N = 1024};

struct Mesh {
    int nt, nv;
    float *r, *c;
    int *t;
};

static float
sq(float x)
{
    return x*x;
}

static float
tri_area(float *a, float *b, float *c)
{
    enum {X, Y, Z};
    double bx, by, bz, cx, cy, cz, A;

    bx = b[X] - a[X];
    by = b[Y] - a[Y];
    bz = b[Z] - a[Z];
    cx = c[X] - a[X];
    cy = c[Y] - a[Y];
    cz = c[Z] - a[Z];
    A = sq(by*cz-bz*cy) + sq(bz*cx-bx*cz) + sq(bx*cy-by*cx);
    return sqrt(A)/2;
}

static float
tri_volume(float *a, float *b, float *c)
{
    enum {X, Y, Z};
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

    V = (ax*by-ay*bx)*cz+(az*bx-ax*bz)*cy+(ay*bz-az*by)*cx;
    return V/6;
}

static int
swap(int n, int size, void *p0) {
    int i;
    char *p, t;
    p = p0;
    while (n--) {
	for (i = 0; i < size/2; i++) {
	    t = p[i];
	    p[i] = p[size - i - 1];
	    p[size - i - 1] = t;
	}
	p += size;
    }
    return 0;
}


static int
eq(const char *a, const char *b)
{
    return strncmp(a, b, N) == 0;
}

static int
line(char *s, FILE *f)
{
    int n;
    if (fgets(s, N, f) == NULL)
	return 1;
    n = strlen(s);
    if (n > 0 && s[n - 1] == '\n')
	s[n - 1] = '\0';
    return 0;
}

static int
scalar(FILE *f, int n, char *name, float *c)
{
    char s[N];
    line(s, f);
    if (sscanf(s, "SCALARS %s float", name) != 1) {
	fprintf(stderr, "%s:%d: expect SCALARS, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);
    }

    line(s, f);
    if (!eq(s, "LOOKUP_TABLE default")) {
	fprintf(stderr, "%s:%d: expect LOOKUP_TABLE, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);
    }
    if (fread(c, sizeof(*c), n, f) != n) {
	fprintf(stderr, "%s:%d: failt to read, n = %d\n", __FILE__, __LINE__, n);
	exit(2);
    }
    return 0;
}

static int
get1(float *r, int i, /**/ float *a)
{
    enum {X, Y, Z};
    a[X] = r[3*i + X]; a[Y] = r[3*i + Y]; a[Z] = r[3*i + Z];
    return 0;
}

static int
get3(struct Mesh *mesh, int t, /**/ float *a, float *b, float *c)
{
  int nt, *tri, i, j, k;
    float *r;
    nt = mesh->nt;
    tri = mesh->t;
    r = mesh->r;
    if (t > nt) {
	fprintf(stderr, "%s:%d: t=%d > nt=%d\n", __FILE__, __LINE__, t, nt);
	exit(1);
    }

    i = tri[3*t];
    j = tri[3*t + 1];
    k = tri[3*t + 2];
    get1(r, i, a);
    get1(r, j, b);
    get1(r, k, c);
    return 0;
}

int
main()
{
    FILE *f = stdin;
    char s[N], name[N];
    int nv, nt;
    float *r, *c;
    int *t, *a, *b;
    int i;
    struct Mesh mesh;
    float u[3], v[3], w[3];

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
	fprintf(stderr, "%s:%d: expect DATASET POLYDATA, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);
    }

    line(s, f);
    if (sscanf(s, "POINTS %d float", &nv) != 1) {
	fprintf(stderr, "%s:%d: expect POINTS, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);
    }
    r = malloc(3*nv*sizeof(*r));
    if (r == NULL) {
	fprintf(stderr, "%s:%d: failt to alloc, nv = %d\n", __FILE__, __LINE__, nv);
	exit(2);
    }
    if (fread(r, sizeof(*r), 3*nv, f) != 3*nv) {
	fprintf(stderr, "%s:%d: failt to read, nv = %d\n", __FILE__, __LINE__, nv);
	exit(2);
    }
    swap(3*nv, sizeof(*r), r);
    while (line(s, f) == 0 && s[0] == '\0') ;
    if (sscanf(s, "POLYGONS %d %*d", &nt) != 1) {
	fprintf(stderr, "%s:%d: expect POLYGONS, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);
    }
    t = malloc(4*nt*sizeof(*t));
    if (t == NULL) {
	fprintf(stderr, "%s:%d: failt to alloc, nt = %d\n", __FILE__, __LINE__, nt);
	exit(2);
    }
    if (fread(t, sizeof(*t), 4*nt, f) != 4*nt) {
	fprintf(stderr, "%s:%d: failt to read, nt = %d\n", __FILE__, __LINE__, nt);
	exit(2);
    }
    for (i = 0, a = b = t; i < nt; i++) {
	b++;
	*a++ = *b++;
	*a++ = *b++;
	*a++ = *b++;
    }
    swap(3*nt, sizeof(*t), t);
    while (line(s, f) == 0 && s[0] == '\0') ;
    if (sscanf(s, "CELL_DATA %*d") != 0) {
	fprintf(stderr, "%s:%d: expect CELL_DATA, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);
    }
    c = malloc(nt*sizeof(*c));
    for (;;) {
      scalar(f, nt, name, c);
      if (eq(name, "cl")) break;
    }
    swap(nt, sizeof(*c), c);

    mesh.r = r;
    mesh.t = t;
    mesh.c = c;
    mesh.nt = nt;
    mesh.nv = nv;

    for (i = 0; i < nt; i++) {
	get3(&mesh, i, u, v, w);
	printf("%ld %.16g %.16g\n",
	       (long)c[i], tri_area(u, v, w), tri_volume(u, v, w));
    }

    free(t);
    free(r);
    free(c);
}
