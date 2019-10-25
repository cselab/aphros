#include <stdio.h>
#include <vtk.h>

int
main(void)
{
  struct VTK *vtk;
  vtk = vtk_read(stdin);
  vtk_fin(vtk);
}
