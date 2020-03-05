#include <stdio.h>
#include <common.h>

int
main(int argc, char **argv)
{
    Point point;
    scalar s;
    coord p;

    point.i = 0;
    point.j = 0;
    point.k = 0;
    point.level = 0;

    p.x = 0;
    p.y = 0;
    p.z = 0;

    init_grid(8);
    printf("%.16g\n", embed_interpolate(point, s, p));
}
