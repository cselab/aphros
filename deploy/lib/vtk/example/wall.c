#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static char me[] = "vtk/wall";
static void usg() {
  fprintf(stderr, "%s < VTK > VTK\n", me);
  exit(1);
}

int main(int argc, char** argv) {
  struct VTK* vtk;
  int n, i, *flag;
  float* cl;

  USED(argc);
  while (*++argv != NULL && *argv[0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(1);
    }

  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }

  n = vtk_nt(vtk);
  flag = malloc(n * sizeof(*flag));
  if (flag == NULL) {
    fprintf(stderr, "%s: malloc failed\n", me);
    exit(2);
  }
  cl = vtk_data(vtk, "cl");
  for (i = 0; i < n; i++)
    flag[i] = (cl[i] == -1 || cl[i] == 0);
  vtk_remove_tri(vtk, flag);
  vtk_write(vtk, stdout);
  free(flag);
  vtk_fin(vtk);
}
