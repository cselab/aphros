#include <stdio.h>
#include <common.h>

int
main(int argc, char **argv)
{
    double c;
    coord n;

    c = 0.25;
    n.x = 0.1;
    n.y = -0.2;
    n.z = -1.8;
    
    printf("%.16g\n", line_alpha(c, n));
}
