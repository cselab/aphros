#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static char me[] = "vtk/field";
int main(int argc, char** argv) {
  USED(argc);
  struct VTK* vtk;
  int nv, nt, n, i, location, type;
  char* name;
  float* df;
  double* dd;
  int* di;
  void* data;

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
  if ((i = vtk_index(vtk, name)) == -1) {
    fprintf(stderr, "%s: no field '%s'\n", me, name);
    exit(2);
  }
  nv = vtk_nv(vtk);
  nt = vtk_nt(vtk);
  data = vtk->data[i];
  location = vtk->location[i];
  type = vtk->type[i];

  n = (location == VTK_CELL) ? nt : nv;

  switch (type) {
    case VTK_FLOAT:
      df = data;
      for (i = 0; i < n; i++)
        printf("%.16g\n", df[i]);
      break;
    case VTK_DOUBLE:
      dd = data;
      for (i = 0; i < n; i++)
        printf("%.16g\n", dd[i]);
      break;
    case VTK_INT:
      di = data;
      for (i = 0; i < n; i++)
        printf("%d\n", di[i]);
      break;
    default:
      fprintf(stderr, "%s: unknown type", me);
      exit(2);
  }
}
