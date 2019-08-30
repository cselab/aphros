#include <stdio.h>
#include <march.h>

enum {X, Y, Z};
static double r = 0.25;
static double lo = -0.5, hi = 0.5;
static int m = 100;

static double
sq(double x)
{
    return x*x;
}

static double
f(double x, double y, double z)
{
    double ans;
    ans = sq(x) + sq(2*y) + sq(3*z) - sq(r);
    return ans;
}

static double O[][3] = {
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {1, 1, 0},
    {0, 0, 1},
    {1, 0, 1},
    {0, 1, 1},
    {1, 1, 1},
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
    fprintf(stderr, "n: %d\n", n);
    for (i = j = 0; i < 3*n; i++) {
	a = d*tri[j++] + x;
	b = d*tri[j++] + y;
	c = d*tri[j++] + z;
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
    static double tri[3*3*5];
    int n, i, j, k, l;
    int u, v, w;
    double x, y, z, d;
    //double cube[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    double cube[8];
    double *o;
    
    printf("# File type: ASCII OBJ\n");

    d = (hi - lo)/(m - 1);
    for (i = 0; i < m; i++)
	for (j = 0; j < m; j++)
	    for (k = 0; k < m; k++) {
		x = lo + d*i;
		y = lo + d*j;
		z = lo + d*k;
		//fprintf(stderr, "%g %g %g\n", x, y, z);
		for (l = 0; l < 8; l++) {
		    o = O[l];
		    cube[l] = f(x + d*o[X], y + d*o[Y], z + d*o[Z]);
		}
		march_cube(cube, &n, tri);
		write(x, y, z, d, n, tri);
	    }
}
