#include <stdio.h>
#include <tgmath.h>
#include "lib.h"

//GetOffset finds the approximate point of intersection of the surface
// between two points with the values Value1 and Value2
float
GetOffset(float a, float b, float x)
{
    double d;
    d = b - a;
    if (d == 0.0)
	return 0.5;
    return (x - a) / d;
}

static int
color0(float x, float y, float z)
{
    return (x > 0.0 ? x : 0.0) + (y < 0.0 ? -0.5 * y : 0.0) + (z < 0.0 ? -0.5 * z : 0.0);
}

//GetColor generates a color from a given position and normal of a point
void
GetColor(struct Vec *Color, struct Vec *n)
{
    Color->x = color0(n->x, n->y, n->z);
    Color->y = color0(n->y, n->z, n->x);
    Color->z = color0(n->z, n->x, n->y);
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

static void
Normalize(struct Vec *v)
{
    double len;
    len = sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
    if (len != 0.0) {
	v->x /= len;
	v->y /= len;
	v->z /= len;
    }
}

//GetNormal() finds the gradient of the scalar field at a point
//This gradient can be used as a very accurate vertx normal for lighting calculations
void
GetNormal(struct Vec *v, float x, float y, float z,
	  double(*f)(double, double, double, void*), void *p)
{
#define F(x, y, z) f((x), (y), (z), p)
    v->x = F(x - 0.01, y, z) - F(x + 0.01, y, z);
    v->y = F(x, y - 0.01, z) - F(x, y + 0.01, z);
    v->z = F(x, y, z - 0.01) - F(x, y, z + 0.01);
    Normalize(v);
#undef F
}
