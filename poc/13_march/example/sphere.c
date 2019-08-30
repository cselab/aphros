#include <stdio.h>
#include <stdlib.h>
#include <march.h>

#define SIZE(x) (sizeof(x)/sizeof(*(x)))

enum { X, Y, Z };
static double r = 0.25;
static double lo = -0.5, hi = 0.5;
static int m = 100;

static double
sq(double x)
{
    return x * x;
}

static double
f(double x, double y, double z)
{
    double ans;

    ans = sq(x) + sq(2*y) + sq(3*z) - sq(r);
    return ans;
}

static double O[][3] = {
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 0, 1, 0 },
    { 1, 1, 0 },
    { 0, 0, 1 },
    { 1, 0, 1 },
    { 0, 1, 1 },
    { 1, 1, 1 },
};

void
    static
write(double x, double y, double z, double d, int n, double *tri)
{
    int i, j, u, v, w;
    double a, b, c;
    static int J = 1;

    if (n == 0)
	return;
    for (i = j = 0; i < 3 * n; i++) {
	a = d * tri[j++] + x;
	b = d * tri[j++] + y;
	c = d * tri[j++] + z;
	printf("v %g %g %g\n", a, b, c);
    }
    for (i = 0; i < n; i++) {
	u = J++;
	v = J++;
	w = J++;
	printf("f %d %d %d\n", u, v, w);
    }
}

int
main(int argc, char **argv)
{
    double tri[3 * 3 * MARCH_NTRI];
    int (*algorithm)(double*, int*, double*);
    int n, i, j, k, l;
    int u, v, w;
    double x, y, z, d;
    double cube[8];
    double *o;
    int stat[MARCH_NTRI] = { 0 };

    if (argv[1] == NULL || argv[1][0] != '-') {
	fprintf(stderr, "%s: needs -c (cube) or -t (tetrahedron)\n", argv[0]);
	exit(2);
    }

    m = 10;
    algorithm = march_cube;
    while (argv[1] != NULL && argv++[1][0] == '-')
	switch (argv[0][1]) {
	case 'c':
	    algorithm = march_cube;
	    break;
	case 't':
	    algorithm = march_tetrahedron;
	    break;
	case 'n':
	    if (argv[1] == NULL) {
		fprintf(stderr, "%s: -n needs an argument\n", argv[0]);
		exit(2);
	    }
	    m = atoi(argv++[1]);
	    break;
	default:
	    fprintf(stderr, "%s: unknow option\n", argv[0]);
	    exit(2);
	}


    printf("# File type: ASCII OBJ\n");
    d = (hi - lo) / (m - 1);
    for (i = 0; i < m; i++)
	for (j = 0; j < m; j++)
	    for (k = 0; k < m; k++) {
		x = lo + d * i;
		y = lo + d * j;
		z = lo + d * k;
		for (l = 0; l < 8; l++) {
		    o = O[l];
		    cube[l] = f(x + d * o[X], y + d * o[Y], z + d * o[Z]);
		}
		algorithm(cube, &n, tri);
		stat[n]++;
		write(x, y, z, d, n, tri);
	    }
    for (i = 0; i < SIZE(stat); i++)
	if (stat[i] > 0)
	    fprintf(stderr, "%2d %7d\n", i, stat[i]);
}
