#include <stdio.h>
#include <stdlib.h>
#include <vtk.h>

static char me[] = "vtk/write";
int
main(void)
{
  struct VTK *vtk;
  int nv, nt, nf, i;
  float *cl;

  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }

  nv = vtk_nv(vtk);
  nt = vtk_nt(vtk);
  nf = vtk_nf(vtk);
  cl = vtk_field(vtk, "cl");
  if (cl == NULL) {
    fprintf(stderr, "%s: cl == NULL\n", me);
    exit(2);
  }

  for (i = 0; i < nt; i++)
    if (cl[i] != -1 && cl[i] != 0)
      printf("%g\n", cl[i]);
  vtk_fin(vtk);
}
