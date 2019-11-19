#include <stdio.h>
#include <stdlib.h>
#include <vtk.h>

#define	USED(x)		if(x);else{}
static char me[] = "vtk/xyz";
static void
usg()
{
  fprintf(stderr, "%s < vtk\n", me);
  exit(1);
}

int
main(int argc, char **argv)
{
  USED(argc);
  struct VTK *vtk;
  int nv, nt, n, i, location, type;
  char *name;
  float *df;
  double *dd;
  int *di;
  void *data;

  while (*++argv != NULL && argv[0][0] == '-')
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
  nv = vtk_nv(vtk);
  for (i = 0; i < nv; i++)
    printf("%.16g %.16g %.16g\n", vtk->x[i], vtk->y[i], vtk->z[i]);

  vtk_fin(vtk);
}
