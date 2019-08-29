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

