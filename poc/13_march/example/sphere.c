#include <stdio.h>
#include <march.h>

#define SIZE(x) (sizeof(x)/sizeof(*(x)))

enum { X, Y, Z };
static double r = 0.25;
static double lo = -0.5, hi = 0.5;
static int m = 20;

static double
sq(double x)
{
    return x * x;
}

static double
f(double x, double y, double z)
{
    double ans;

    ans = sq(x) + sq(2 * y) + sq(3 * z) - sq(r);
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
main()
{
    static double tri[3*3*MARCH_NTRI];
    int n, i, j, k, l;
    int u, v, w;
    double x, y, z, d;
    double cube[8];
    double *o;
    int stat[MARCH_NTRI] = { 0 };
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
		march_tetrahedron(cube, &n, tri);
		stat[n]++;
		write(x, y, z, d, n, tri);
	    }
    for (i = 0; i < SIZE(stat); i++)
	if (stat[i] > 0)
	    fprintf(stderr, "%2d %7d\n", i, stat[i]);
}
