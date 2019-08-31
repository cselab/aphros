#include <stdio.h>
#include "mc.h"

static int nx = 60, ny = 60, nz = 60;

void compute_data(MarchingCubes & mc);
int
main()
{
    MarchingCubes mc(nx, ny, nz);
    compute_data(mc);
    mc.run();
    mc.clean_temps();
    writeObj(&mc);
    mc.clean_all();
    return 0;
}

void
compute_data(MarchingCubes & mc)
{
    double x, y, z;
    double sx, sy, sz;
    double tx, ty, tz;
    double r;

    r = 1.85;
    sx = (double) nx / 16;
    sy = (double) ny / 16;
    sz = (double) nz / 16;
    tx = (double) nx / (2 * sx);
    ty = (double) ny / (2 * sy) + 1.5f;
    tz = (double) nz / (2 * sz);
    for (int k = 0; k < nz; k++) {
	z = ((double) k) / sz - tz;
	for (int j = 0; j < ny; j++) {
	    y = ((double) j) / sy - ty;
	    for (int i = 0; i < nx; i++) {
		x = ((double) i) / sx - tx;
		mc.set_data((double) (x*x + y*y + z*z - r*r), i, j, k);
	    }
	}
    }
}
