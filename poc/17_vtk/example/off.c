#include <stdio.h>
#include <stdlib.h>
#include <vtk.h>
#define	USED(x)		if(x);else{}
static char me[] = "vtk/off";
static void
usg()
{
    fprintf(stderr, "%s < VTK > OFF\n", me);
    exit(1);
}

int
main(int argc, char **argv)
{
  struct VTK *vtk;

  USED(argc);
  while (*++argv != NULL && argv[0][0] == '-' )
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
  vtk_off(vtk, stdout);
  vtk_fin(vtk);
}
