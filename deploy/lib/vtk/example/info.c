#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

static char me[] = "vtk/info";
int main(void) {
  struct VTK* vtk;
  int nv, nt, nf, i, location, type;

  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }
  nv = vtk_nv(vtk);
  nt = vtk_nt(vtk);
  nf = vtk_nf(vtk);

  printf("nv %d\n", nv);
  printf("nt %d\n", nt);
  printf("nf %d\n", nf);
  for (i = 0; i < nf; i++) {
    location = vtk->location[i];
    type = vtk->type[i];
    printf(
        "%s %s %s %d\n", vtk->name[i],
        type == VTK_FLOAT    ? "float"
        : type == VTK_DOUBLE ? "double"
        : type == VTK_INT    ? "int"
                             : "unknown",
        location == VTK_CELL ? "cell" : "point", vtk->rank[i]);
  }
  vtk_fin(vtk);
}
