#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {N = 1024};

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
    int nv, nt, nt0;
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
    if (fread(r, sizeof(*r), 3*nv, f) != 3*nv) {
	fprintf(stderr, "%s:%d: failt to read, nv = %d\n", __FILE__, __LINE__, nv);
	exit(2);
    }
    swap(3*nv, sizeof(*r), r);
    //fprintf(stderr, "%g %g %g\n", r[0], r[1], r[2]);
    //fprintf(stderr, "%g %g %g\n", r[3*nv-3], r[3*nv-2], r[3*nv-1]);

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
    fprintf(stderr, "%d %d %d\n", t[0], t[1], t[2]);
    fprintf(stderr, "%d %d %d\n", t[3*nt-3], t[3*nt-2], t[3*nt-1]);
    while (line(s, f) == 0 && s[0] == '\0') ;
    if (sscanf(s, "CELL_DATA %*d") != 0) {
	fprintf(stderr, "%s:%d: expect CELL_DATA, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);
    }
    line(s, f);
    if (!eq(s, "SCALARS c float")) {
	fprintf(stderr, "%s:%d: expect SCALARS c, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);	
    }
    line(s, f);
    if (!eq(s, "LOOKUP_TABLE default")) {
	fprintf(stderr, "%s:%d: expect LOOKUP_TABLE, got '%s'\n", __FILE__, __LINE__, s);
	exit(2);	
    }    
    c = malloc(nt*sizeof(*c));
    if (c == NULL) {
	fprintf(stderr, "%s:%d: failt to alloc, nt = %d\n", __FILE__, __LINE__, nt);
	exit(2);
    }
    if (fread(c, sizeof(*c), nt, f) != nt) {
	fprintf(stderr, "%s:%d: failt to read, nt = %d\n", __FILE__, __LINE__, nt);
	exit(2);
    }
    swap(nt, sizeof(*c), c);
    for (i = 0; i < nt; i++)
	fprintf(stdout, "%f\n", c[i]);
    free(t);
    free(r);
    free(c);
}
