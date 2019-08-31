#include <stdlib.h>
#include <stdio.h>
#include "mc.h"

static int nx = 60, ny = 60, nz = 60;

int
main()
{
    int i, j, k, m;
    double x, y, z;
    double sx, sy, sz;
    double tx, ty, tz;
    double r;
    double *buf;
    MarchingCubes mc(nx, ny, nz);
    buf = (double*)malloc(nx*ny*nz*sizeof(*buf));
    if (buf == NULL) {
	fprintf(stderr, "alloc failed\n");
	exit(2);
    }

    r = 1.85;
    sx = nx / 16.0;
    sy = ny / 16.0;
    sz = nz / 16.0;
    tx = nx / (2 * sx);
    ty = ny / (2 * sy);
    tz = nz / (2 * sz);
    for (m = k = 0; k < nz; k++)
	for (j = 0; j < ny; j++)
	    for (i = 0; i < nx; i++) {
		x = i / sx - tx;
		y = j / sy - ty;
		z = k / sz - tz;
		buf[m++] = x * x + y * y + z * z - r * r;
		mc.set_data(x * x + y * y + z * z - r * r, i, j, k);
	    }

    mc.run();
    mc.clean_temps();
    writeObj(&mc);
    mc.clean_all();
    free(buf);
    return 0;
}

