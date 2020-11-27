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
  int i;
  int Invert;
  int nt;
  int nv;
  int Pflag;
  int tmp;
  int* tri;
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
  printf("size = %g\n", info.size);
  printf("nx = %d\n", info.nx);
  printf("ny = %d\n", info.ny);  
  inside_fin(inside);
  inside_mesh_fin(tri, ver);
}
