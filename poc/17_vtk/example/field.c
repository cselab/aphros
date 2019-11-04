#include <stdio.h>
#include <stdlib.h>
#include <vtk.h>

#define	USED(x)		if(x);else{}
static char me[] = "vtk/field";
int
main(int argc, char **argv)
{
  USED(argc);
  struct VTK *vtk;
  int nv, nt, n, i, location;
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
  if ((i = vtk_index(vtk, name)) != -1) {
    fprintf(stderr, "%s: no field '%s'\n", me, name);
    exit(2);
  }
  nv = vtk_nv(vtk);
  nt = vtk_nt(vtk);
  data = vtk->data[i];
  location = vtk->location[i];

  n = (location == VTK_CELL) ? nt : nv;
  for (i = 0; i < n; i++) {
    printf("%g\n", data[i]);
  }
}
