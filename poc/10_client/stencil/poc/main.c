#include <stdlib.h>
#include <stdio.h>
#include "lib.h"

static int a[10*10*10];

static int
vtk(const char *path, int n, const int *a)
{
    FILE *f;
    int i;
    double spa, org;
    const char *s = \
	"# vtk DataFile Version 2.0\n"	\
	"generated with poc\n"		\
	"ASCII\n"			\
	"DATASET STRUCTURED_POINTS\n"	\
	"DIMENSIONS %d %d %d\n"		\
	"ORIGIN %.16g %.16g %.16g\n"	\
	"SPACING %.16g %.16g %.16g\n"	\
	"CELL_DATA %d\n"		\
	"SCALARS color int\n"		\
	"LOOKUP_TABLE default\n"	\
	;

    if (n <= 1) {
	fprintf(stderr, "%s:%d: n=%d <= 1\n", __FILE__, __LINE__, n);
	return 1;
    }
    spa = 1.0/(n - 1);
    org = 0.5 - spa/2;
    
    f = fopen(path, "w");
    if (!f) {
	fprintf(stderr, "%s:%d: can't open '%s'\n", __FILE__, __LINE__, path);
	return 1;
    }
    fprintf(f, s, n + 1, n + 1, n + 1, org, org, org, spa, spa, spa, n*n*n);
    for (i = 0; i < n*n*n; i++)
	fprintf(f, "%d\n", a[i]);
    fclose(f);
    return 0;
}

int
main()
{
    
    const char *in = "a.vtk", *out = "b.vtk";
    int i, n, cnt, ret;

    ret = scanf("%d", &n);
    if (ret != 1) {
	fprintf(stderr, "can't read n\n");
	exit(2);
    }
	
    
    for (i = 0; i < n*n*n; i++) {
	ret = scanf("%d", &a[i]);
	if (ret != 1) {
	    fprintf(stderr, "can't read a[%d]\n", i);
	    exit(2);
	}
    }
    
    vtk(in, n, a);
    color(n, &cnt, a);
    vtk(out, n, a);

    fprintf(stderr, "%s\n", in);
    fprintf(stderr, "%s\n", out);
    fprintf(stderr, "cnt: %d\n", cnt);
}
