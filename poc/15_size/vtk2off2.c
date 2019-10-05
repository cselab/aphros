#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include "arg.h"

enum { N = 1024 };
static char me[] = "ch.vtk2off";

char *argv0;
static int nv, nt;
static float *r;
static int *t, *c;
static float *cl;

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

#define FWRITE(n, p, f) \
    do {								\
	if ((int)fwrite(p, sizeof(*(p)), (n), (f)) != (n)) {		\
	    fprintf(stderr, "%s:%d: failt to write, n = %d\n", __FILE__, __LINE__, n); \
	    exit(2);							\
	}								\
} while(0)
#define SWAP(n, p) swap(n, sizeof(*(p)), p)
#define HASH(x, s) hash(x, sizeof(*(x)), s)

static void
usg(void)
{
    fprintf(stderr, "usage: %s < VTK > OFF\n", me);
    exit(0);
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
read_vtk(void)
{
    FILE *f;
    char s[N], name[N];
    int i, *a, *b, *t0;

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
    SWAP(3 * nv, r);

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
    SWAP(3*nt, t);
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
    SWAP(nt, cl);
    return 0;
}

static int
write_off()
{
    FILE *f;
    int i, j, k, *t0;
    float *r0;
    f = stdout;
    fputs("OFF BINARY\n", f);
    MALLOC(5*nt, &t0);
    MALLOC(3*nv, &r0);
    for (i = 0; i < 3*nv; i++)
	r0[i] = r[i];
    SWAP(3*nv, r0);
    for (i = j = k = 0; i < nt; i++) {
	t0[j++] = 3;
	t0[j++] = t[k++];
	t0[j++] = t[k++];
	t0[j++] = t[k++];
	t0[j++] = 0;
    }
    SWAP(5*nt, t0);
    i = nv; SWAP(1, &i); FWRITE(1, &i, f);
    i = nt; SWAP(1, &i); FWRITE(1, &i, f);
    i =  0; SWAP(1, &i); FWRITE(1, &i, f);
    FWRITE(3*nv, r0, f);
    FWRITE(5*nt, t0, f);
    free(t0);
    free(r0);
    return 0;
}

static int
wall(void)
{
    int i, j;
    for (i = j = 0; i < nt; i++)
    {
	if ((int)cl[i] != -1) {
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
write_vtk(void)
{
    FILE *f;
    int i, j, k, *t0;
    float *r0;

    f = stdout;
    MALLOC(4*nt, &t0);
    MALLOC(3*nv, &r0);
    for (i = 0; i < 3*nv; i++)
	r0[i] = r[i];
    for (i = j = k = 0; i < nt; i++) {
	t0[j++] = 3;
	t0[j++] = t[k++];
	t0[j++] = t[k++];
	t0[j++] = t[k++];
    }
    SWAP(3*nv, r0);
    SWAP(4*nt, t0);
    fprintf(f, "# vtk DataFile Version 2.0\n"
	    "Interface from marching cubes\n"
	    "BINARY\n"
	    "DATASET POLYDATA\n"
	    "POINTS %d float\n", nv);
    FWRITE(3*nv, r0, f);
    fprintf(f, "POLYGONS %d %d\n", nt, 4*nt);
    FWRITE(4*nt, t0, f);
    fprintf(f, "CELL_DATA %d\n"
	    "SCALARS c int\n"
	    "LOOKUP_TABLE default\n", nt);
    for (i = 0; i < nt; i++) {
	t0[i] = c[i];
    }
    SWAP(nt, t0);
    FWRITE(nt, t0, f);
    free(t0);
    free(r0);
    return 0;
}

static void
hash(const void *p0, int n, char *s)
{
    const char *p;
    int i;
    p = p0;
    for (i = 0; i < n; i++)
	s[i] = p[i];
    s[n] = '\0';
}
static int
color(int *pnc)
{
    ENTRY e, *p;
    char s[42];
    int i;
    size_t j;
    float c0;

    hcreate(nt);
    MALLOC(nt, &c);
    for (i = j = 0; i < nt; i++) {
	c0 = cl[i];
	HASH(&c0, s);
	e.key = s;
	if (hsearch(e, FIND) == NULL) {
	    e.data = (void*)(j++);
	    hsearch(e, ENTER);
	}
	p = hsearch(e, FIND);
	c[i] = (size_t)(p->data);
    }
    hdestroy();
    *pnc = j;
    return 0;
}

int
main(int argc, char **argv)
{
    int (*Write)(void);
    int nc;

    Write = write_off;
    ARGBEGIN {
	case 'v':
	    Write = write_vtk;
	    break;
	case 'h':
	    usg();
    } ARGEND;

    read_vtk();
    wall();
    color(&nc);
    fprintf(stderr, "nc = %d\n", nc);
    Write();

    free(r);
    free(t);
    free(cl);
}
