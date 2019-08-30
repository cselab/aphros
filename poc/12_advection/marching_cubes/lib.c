#include <stdio.h>
#include <tgmath.h>
#include "lib.h"



double
GetOffset(double a, double b)
{
    double d;

    d = a - b;
    if (d == 0.0)
	return 0.5;
    return a/d;
}

static int
color0(double x, double y, double z)
{
    return (x > 0.0 ? x : 0.0) + (y < 0.0 ? -0.5 * y : 0.0) + (z <
							       0.0 ? -0.5 *
							       z : 0.0);
}
void
GetColor(double x, double y, double z, double *u, double *v, double *w)
{
    *u = color0(x, y, z);
    *v = color0(y, z, x);
    *w = color0(z, x, y);
}

void
PrintHelp(void)
{
    printf
	("Marching Cubes Example by Cory Bloyd (dejaspaminacan@my-deja.com)\n\n");

    printf("+/-  increase/decrease sample density\n");
    printf("PageUp/PageDown  increase/decrease surface value\n");
    printf("s  change sample function\n");
    printf("c  toggle marching cubes / marching tetrahedrons\n");
    printf("w  wireframe on/off\n");
    printf("l  toggle lighting / color-by-normal\n");
    printf("Home  spin scene on/off\n");
    printf("End  source point animation on/off\n");
}

static double
sq(double x)
{
    return x*x;
}

static void
Normalize(double *u, double *v, double *w)
{
    double len;

    len = sqrt(sq(*u) + sq(*v) + sq(*w));
    if (len != 0.0) {
	*u /= len;
	*v /= len;
	*w /= len;
    }
}


void
GetNormal(double x, double y, double z, double *u, double *v, double *w,
	  double (*f) (double, double, double, void *), void *p)
{
#define F(x, y, z) f((x), (y), (z), p)
    *u = F(x - 0.01, y, z) - F(x + 0.01, y, z);
    *v = F(x, y - 0.01, z) - F(x, y + 0.01, z);
    *w = F(x, y, z - 0.01) - F(x, y, z + 0.01);
    Normalize(u, v, w);
#undef F
}
