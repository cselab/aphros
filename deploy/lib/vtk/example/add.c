#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

static char me[] = "vtk/add";
int main(void) {
  struct VTK* vtk;
  int n, i;
  double* id;

  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }

  n = vtk_nt(vtk);
  vtk_add(vtk, "id", VTK_CELL, VTK_DOUBLE);
  id = vtk_data(vtk, "id");
  for (i = 0; i < n; i++)
    id[i] = i;
  vtk_write(vtk, stdout);
  vtk_fin(vtk);
}
