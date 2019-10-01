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
    int nv;
    float *r;

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
    fprintf(stderr, "%g %g %g\n", r[0], r[1], r[2]);
    fprintf(stderr, "%g %g %g\n", r[3*nv-3], r[3*nv-2], r[3*nv-1]);

    
    free(r);
}
