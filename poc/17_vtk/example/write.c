#include <stdio.h>
#include <stdlib.h>
#include <vtk.h>

static char me[] = "vtk/write";
int
main(void)
{
  struct VTK *vtk;
  int nv, nt, i;

  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }

  nv = vtk_nv(vtk);
  nt = vtk_nt(vtk);

  vtk_fin(vtk);
}
