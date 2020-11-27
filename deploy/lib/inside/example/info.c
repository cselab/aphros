/*
c99 info.c -I.. -L.. -linside -lm
*/
#include <stdio.h>
#include <stdlib.h>
#include <inside.h>

const char* me = "info";
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static void usg(void) {
  fprintf(stderr, "%s < off\n", me);
  exit(1);
}

int main(int argc, const char** argv) {
  USED(argc);
  enum { X, Y, Z };
  double* ver;
  int nt;
  int nv;
  int* tri;
  double lo[3];
  double hi[3];
  struct Inside* inside;
  struct InsideInfo info;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(2);
    }
  if (argv[0] == NULL) {
    fprintf(stderr, "%s: mesh file is missing\n", me);
    exit(2);
  }
  if (inside_mesh_read(argv[0], &nt, &tri, &nv, &ver) != 0) {
    fprintf(stderr, "%s: fail to read mesh '%s'\n", me, argv[0]);
    exit(2);
  }
  inside_ini(nt, tri, ver, &inside);
  inside_info(inside, &info);
  inside_box(inside, lo, hi);
  printf("box: [%g %g %g] [%g %g %g]\n",
	 lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
  printf("size = %g\n", info.size);
  printf("nx = %d\n", info.nx);
  printf("ny = %d\n", info.ny);  
  inside_fin(inside);
  inside_mesh_fin(tri, ver);
}
