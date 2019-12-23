#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

static char me[] = "vtk/write";
int main(void) {
  struct VTK* vtk;

  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }
  vtk_write(vtk, stdout);
  vtk_fin(vtk);
}
