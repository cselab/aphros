#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {N = 1024};

static int *u_root;
static void
u_ini(int n)
{
    u_root = malloc(n*sizeof(*u_root));
    if (u_root == NULL) {
	fprintf(stderr, "%s:%d: alloc failed, n = %d\n", __FILE__, __LINE__, n);
	exit(2);
    }
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

int
main()
{
    FILE *f = stdin;
    char s[N];
    int nv, nt, u, v, w;
    float *r, *c;
    int *t, *a, *b;
    int i;

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
    if ((int)fread(r, sizeof(*r), 3*nv, f) != 3*nv) {
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
    if ((int)fread(t, sizeof(*t), 4*nt, f) != 4*nt) {
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
    u_ini(nv);
    for (i = 0 ; i < nt; i++) {
	u = t[3*i];
	v = t[3*i + 1];
	w = t[3*i + 2];
	u_union(u, v);
	u_union(v, w);
    }

    c = malloc(nt*sizeof(*c));
    if (c == NULL) {
	fprintf(stderr, "%s:%d: failt to alloc, nt = %d\n", __FILE__, __LINE__, nt);
	exit(2);
    }
    for (i = 0; i < nt; i++) {
	u = t[3*i];
	c[i] = u_find(i);
    }

    
    
    u_fin();
    free(c);
    free(t);
    free(r);
}
