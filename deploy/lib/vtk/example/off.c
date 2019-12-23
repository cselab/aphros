#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static char me[] = "vtk/off";
static void usg() {
  fprintf(stderr, "%s [-f field] < VTK > OFF\n", me);
  exit(1);
}

int main(int argc, char** argv) {
  struct VTK* vtk;
  const char* field;

  USED(argc);
  field = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'f':
        argv++;
        if ((field = *argv) == NULL) {
          fprintf(stderr, "%s: -f needs an argument\n", me);
          exit(2);
        }
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

  if (field == NULL)
    vtk_off(vtk, stdout);
  else
    vtk_off_color(vtk, field, stdout);
  vtk_fin(vtk);
}
