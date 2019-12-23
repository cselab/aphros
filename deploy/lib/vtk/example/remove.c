#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static char me[] = "vtk/remove";

static void usg() {
  fprintf(stderr, "%s [field ..] < VTK > VTK\n", me);
  exit(1);
}

int main(int argc, char** argv) {
  USED(argc);
  struct VTK* vtk;
  int status;

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

  while (*argv != NULL) {
    status = vtk_remove(vtk, *argv);
    if (status != 0) fprintf(stderr, "%s: no field '%s'\n", me, *argv);
    argv++;
  }

  vtk_write(vtk, stdout);
  vtk_fin(vtk);
}
