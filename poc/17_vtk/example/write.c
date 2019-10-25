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

  for (i = 0; i < nt; i++)
      fprintf(stderr, "%d %d %d\n", vtk->t0[i], vtk->t1[i], vtk->t2[i]);

  

  vtk_fin(vtk);
}
