#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static char me[] = "vtk/one";

static void usg() {
  fprintf(stderr, "%s -f field -k [int] VTK > VTK\n", me);
  fprintf(stderr, "extract one bubble\n");
  exit(1);
}

int main(int argc, char** argv) {
  struct VTK* vtk;
  const char* name;
  int nt, i, k, location, type, n;
  float* df;
  int *di, Key, Field;
  double* dd;
  void* data;

  USED(argc);
  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }
  Key = Field = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'f':
        argv++;
        if ((name = argv[0]) == NULL) {
          fprintf(stderr, "%s: -f needs an argument\n", me);
          exit(2);
        }
        if ((i = vtk_index(vtk, name)) == -1) {
          fprintf(stderr, "%s: no field '%s'\n", me, name);
          exit(2);
        }
        if (vtk->location[i] != VTK_CELL) {
          fprintf(stderr, "%s: %s is not a cell data\n", me, name);
          exit(1);
        }
        if (vtk->rank[i] != VTK_SCALAR) {
          fprintf(stderr, "%s: %s is not a scalar\n", me, name);
          exit(1);
        }
        Field = 1;
        break;
      case 'k':
        argv++;
        if ((name = argv[0]) == NULL) {
          fprintf(stderr, "%s: -f needs an argument\n", me);
          exit(2);
        }
        k = atoi(name);
        Key = 1;
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(1);
    }
  if (!Field) {
    fprintf(stderr, "%s: field (-f) is not given\n", me);
    exit(1);
  }
  if (!Key) {
    fprintf(stderr, "%s: key (-k) is not given\n", me);
    exit(1);
  }
  data = vtk->data[i];
  type = vtk->type[i];
  nt = vtk_nt(vtk);
  switch (type) {
    case VTK_FLOAT:
      df = data;
      break;
    case VTK_DOUBLE:
      dd = data;
      break;
    case VTK_INT:
      di = data;
      break;
    default:
      fprintf(stderr, "%s: unknown type", me);
      exit(2);
  }
}
