#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h>
#include "mc.h"

#define T2 ( \
	((x*x+y*y+z*z+R*R-r*r)*(x*x+y*y+z*z+R*R-r*r)-4*R*R*(x*x+y*y))*	\
	((x*x+(y+R)*(y+R)+z*z+R*R-r*r)*(x*x+(y+R)*(y+R)+z*z+R*R-r*r)-4*R*R*((y+R)*(y+R)+z*z)) \
	)
#define MC (\
    - 26.5298*(1-x)*(1-y)*(1-z) + 81.9199*x*(1-y)*(1-z) - 100.68*x*y*(1-z) + 3.5498*(1-x)*y*(1-z) \
    + 24.1201*(1-x)*(1-y)*  z   - 74.4702*x*(1-y)*  z   + 91.5298*x*y*  z  - 3.22998*(1-x)*y*  z  \
	)
#define CHAIR (\
	x*x+y*y+z*z-0.95f*25)*(x*x+y*y+z*z-0.95f*25)-0.8f*((z-5)*(z-5)-2*x*x)*((z+5)*(z+5)-2*y*y \
	    )

#define CUSHIN (							\
	z*z*x*x - z*z*z*z - 2*z*x*x + 2*z*z*z + x*x - z*z - (x*x - z)*(x*x - z) - y*y*y*y - 2*x*x*y*y - y*y*z*z + 2*y*y*z + y*y \
	)


int
main()
{
    int nx = 40, ny = 40, nz = 40;
    int i, j, k, m;
    double x, y, z;
    double sx, sy, sz;
    double tx, ty, tz;
    double r, R;
    double *buf;

    buf = (double *) malloc(nx * ny * nz * sizeof(*buf));
    if (buf == NULL) {
	fprintf(stderr, "alloc failed\n");
	exit(2);
    }
    r = 1.85;
    R = 4;
    sx = nx / 16.0;
    sy = ny / 16.0;
    sz = nz / 16.0;
    tx = nx / (2 * sx);
    ty = ny / (2 * sy) + 1.5;
    tz = nz / (2 * sz);
    for (m = k = 0; k < nz; k++)
	for (j = 0; j < ny; j++)
	    for (i = 0; i < nx; i++) {
		x = i / sx - tx;
		y = j / sy - ty;
		z = k / sz - tz;
		buf[m++] = T2;
	    }

    MarchingCubes(nx, ny, nz, buf);
    run();
    writeObj();
    free(buf);
    return 0;
}
