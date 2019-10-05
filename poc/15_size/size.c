#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include "arg.h"

enum { N = 1024 };
static double pi = 3.141592653589793;
static char me[] = "ch.size";
static int UseCl, WriteVolume;

char *argv0;
static int nv, nt;
static float *r;
static int *t, *t0, *c;
static float *cl;
static double *dot;
static double *cx, *cy, *cz, *volume;
static int *cnt;

static void
usg(void)
{
    fprintf(stderr, "usage: %s [-c]\n", me);
    exit(0);
}

#define MALLOC(n, p)							\
    do {								\
	*(p) = malloc(n*sizeof(**(p)));					\
	if (*(p) == NULL) {						\
	    fprintf(stderr, "%s:%d: alloc failed, n = %d\n", __FILE__, __LINE__, n); \
	    exit(2);							\
	}								\
    } while(0)

#define FREAD(n, p, f) \
    do {								\
	if ((int)fread(p, sizeof(*(p)), (n), (f)) != (n)) {		\
	    fprintf(stderr, "%s:%d: failt to read, n = %d\n", __FILE__, __LINE__, n); \
	    exit(2);							\
	}								\
} while(0)

struct Mesh {
    int nt, nv;
    float *r;
    int *t;
};

double
tri_volume_y(float *a, float *b, float *c)
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

    V = (cy+by+ay)*((bz-az)*(cx-ax)-(bx-ax)*(cz-az));
    return V/6;
}


static int
eq(const char *a, const char *b)
{
    return strncmp(a, b, N) == 0;
}

static int
line(char *s, FILE * f)
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
tri_center(float *a, float *b, float *c, /**/ double *x, double *y,
	   double *z)
{
    enum { X, Y, Z };

    *x += (a[X] + b[X] + c[X]) / 3;
    *y += (a[Y] + b[Y] + c[Y]) / 3;
    *z += (a[Z] + b[Z] + c[Z]) / 3;
    return 0;
}

static double
tri_dot(float *a, float *b, float *c, double x, double y, double z)
{
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

static int
scalar(FILE * f, int n, char *name, float *c)
{
    char s[N];

    line(s, f);
    if (sscanf(s, "SCALARS %s float", name) != 1) {
	fprintf(stderr, "%s:%d: expect SCALARS, got '%s'\n", __FILE__,
		__LINE__, s);
	exit(2);
    }
    line(s, f);
    if (!eq(s, "LOOKUP_TABLE default")) {
	fprintf(stderr, "%s:%d: expect LOOKUP_TABLE, got '%s'\n", __FILE__,
		__LINE__, s);
	exit(2);
    }
    FREAD(n, c, f);
    return 0;
}

static int
get1(float *r, int i, /**/ float *a)
{
    enum { X, Y, Z };

    a[X] = r[3 * i + X];
    a[Y] = r[3 * i + Y];
    a[Z] = r[3 * i + Z];
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
	fprintf(stderr, "%s:%d: t=%d > nt=%d\n", __FILE__, __LINE__, t,
		nt);
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

static int *u_root;
static int
max_arg(int n, int *a)
{
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

static void
u_ini(int n)
{
    MALLOC(n, &u_root);
    while (n--)
	u_root[n] = n;
}

static void
u_fin()
{
    free(u_root);
}

static int
u_find(int v)
{
    if (v == u_root[v])
	return v;
    return u_root[v] = u_find(u_root[v]);
}

static void
u_union(int a, int b)
{
    a = u_find(a);
    b = u_find(b);
    if (a != b)
	u_root[b] = a;
}

static int
swap(int n, int size, void *p0)
{
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

static int
read_vtk(void)
{
    FILE *f;
    char s[N], name[N];
    int i;
    int *a, *b;

    f = stdin;
    if (line(s, f) != 0) {
	fprintf(stderr, "%s:%d: failt to read\n", __FILE__, __LINE__);
	exit(2);
    }

    if (!eq(s, "# vtk DataFile Version 2.0")) {
	fprintf(stderr, "%s:%d: not a vtk file: '%s'\n", __FILE__,
		__LINE__, s);
	exit(2);
    }

    line(s, f);
    line(s, f);
    if (!eq(s, "BINARY")) {
	fprintf(stderr, "%s:%d: expect BINARY, got '%s'\n", __FILE__,
		__LINE__, s);
	exit(2);
    }

    line(s, f);
    if (!eq(s, "DATASET POLYDATA")) {
	fprintf(stderr, "%s:%d: expect DATASET POLYDATA, got '%s'\n",
		__FILE__, __LINE__, s);
	exit(2);
    }

    line(s, f);
    if (sscanf(s, "POINTS %d float", &nv) != 1) {
	fprintf(stderr, "%s:%d: expect POINTS, got '%s'\n", __FILE__,
		__LINE__, s);
	exit(2);
    }
    MALLOC(3 * nv, &r);
    FREAD(3 * nv, r, f);
    swap(3 * nv, sizeof(*r), r);

    while (line(s, f) == 0 && s[0] == '\0');
    if (sscanf(s, "POLYGONS %d %*d", &nt) != 1) {
	fprintf(stderr, "%s:%d: expect POLYGONS, got '%s'\n", __FILE__,
		__LINE__, s);
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
    swap(3 * nt, sizeof(*t), t);
    while (line(s, f) == 0 && s[0] == '\0');
    if (sscanf(s, "CELL_DATA %*d") != 0) {
	fprintf(stderr, "%s:%d: expect CELL_DATA, got '%s'\n", __FILE__,
		__LINE__, s);
	exit(2);
    }
    MALLOC(nt, &cl);
    for (;;) {
	scalar(f, nt, name, cl);
	if (eq(name, "cl"))
	    break;
    }
    swap(nt, sizeof(*cl), cl);
    return 0;
}

static int
wall(void)
{
    int i, j;
    for (i = j = 0; i < nt; i++)
    {
	if (cl[i] != -1) {
	    t[3*j] = t[3*i];
	    t[3*j+1] = t[3*i+1];
	    t[3*j+2] = t[3*i+2];
	    cl[j] = cl[i];
	    j++;
	}
    }
    nt = j;
    return 0;
}

static int
color(int *pnb)
{
    int nb, u, v, w;
    int *id;
    int i, j, k;
    struct Mesh mesh;
    float x[3], y[3], z[3];

    MALLOC(nt, &c);
    MALLOC(nt, &id);

    if (UseCl) {
	ENTRY e, *p;
	char s[42];
	size_t j;
	hcreate(nt);
	for (i = j = 0; i < nt; i++) {
	    sprintf(s, "%ld", (size_t)cl[i]);
	    e.key = s;
	    if (hsearch(e, FIND) == NULL) {
		e.data = (void*)(j++);
		hsearch(e, ENTER);
	    }
	    p = hsearch(e, FIND);
	    c[i] = (size_t)(p->data);
	}
	hdestroy();
    } else {
	u_ini(nv);
	for (i = 0; i < nt; i++) {
	    u = t[3 * i];
	    v = t[3 * i + 1];
	    w = t[3 * i + 2];
	    u_union(u, v);
	    u_union(v, w);
	}
	for (i = 0; i < nt; i++) {
	    u = t[3 * i];
	    c[i] = u_find(u);
	}
	u_fin();
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

    k = max_arg(nb, cnt);       /* water = 0 */
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
	if (k != 0 && dot[k] > 0)
	    c[i] = 0;
    }
    free(id);
    *pnb = nb;
    return 0;
}

int
main(int argc, char **argv)
{
    double R, V;
    int i, k, nb;
    struct Mesh mesh;
    float x[3], y[3], z[3];

    UseCl = WriteVolume = 0;
    ARGBEGIN {
	case 'c':
	    UseCl = 1;
	    break;
	case 'v':
	    WriteVolume = 1;
	    break;
	case 'h':
	    usg();
    } ARGEND;

    read_vtk();
    wall();
    color(&nb);
    MALLOC(nb, &volume);
    mesh.r = r;
    mesh.t = t;
    mesh.nt = nt;
    mesh.nv = nv;
    for (i = 0; i < nb; i++)
	volume[i] = 0;
    for (i = 0; i < nt; i++) {
	if (c[i] == -1)
	    continue;
	k = c[i];
	get3(&mesh, i, x, y, z);
	volume[k] += tri_volume_y(x, y, z);
    }

    if (WriteVolume) {
	V = 0;
	for (i = 0; i < nb; i++)
	    if (volume[i] > 0 && dot[i] < 0)
		V += fabs(volume[i]);
	printf("%.16g\n", V);
    } else {
	printf("x y z r n\n");
	for (i = 0; i < nb; i++)
	    if (volume[i] > 0 && dot[i] < 0) {
		V = volume[i];
		R = pow(3 * V / (4 * pi), 1.0 / 3.0);
		printf("%.16e %.16e %.16e %.16e %d\n",
		       cx[i], cy[i], cz[i], R, cnt[i]);
	    }
    }

    free(c);
    free(cl);
    free(cnt);
    free(cx);
    free(cy);
    free(cz);
    free(dot);
    free(r);
    free(t);
    free(t0);
    free(volume);
}
