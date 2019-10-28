#include <stdio.h>
#include <stdlib.h>
#include <vtk.h>

static char me[] = "vtk/field";
int
main(int argc, char **argv)
{
  struct VTK *vtk;
  int nv, nt, n, i, location, type, rank;
  char *name;
  float *data;

  argv++;
  if (argv[0] == NULL) {
    fprintf(stderr, "%s: needs an argument\n", me);
    exit(2);
  }
  name = argv++[0];

  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }
  ;
  if (vtk_index(vtk, name, &i) != 0) {
    fprintf(stderr, "%s: no field '%s'\n", me, name);
    exit(2);
  }
  nv = vtk_nv(vtk);
  nt = vtk_nt(vtk);
  data = vtk->data[i];
  location = vtk->location[i];
  type = vtk->type[i];
  rank = vtk->rank[i];

  n = (location == VTK_CELL) ? nt : nv;
  for (i = 0; i < n; i++) {
    printf("%g\n", data[i]);
  }
}
