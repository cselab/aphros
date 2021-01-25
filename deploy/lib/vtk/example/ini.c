#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

enum { nv = 3, nt = 1 };

// static char me[] = "vtk/ini";
int main(void) {
  int i;
  double* c;
  double x[] = {0, 1, 0}, y[] = {0, 0, 1}, z[] = {0, 0, 0};
  int t0[] = {0}, t1[] = {1}, t2[] = {2};
  struct VTK* vtk;

  vtk = vtk_ini(nv, x, y, z, nt, t0, t1, t2);

  vtk_add(vtk, "c", VTK_CELL, VTK_DOUBLE);
  c = vtk_data(vtk, "c");

  for (i = 0; i < nt; i++)
    c[i] = i;

  vtk_write(vtk, stdout);
  vtk_fin(vtk);
}
