#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inside.h>

static const char* me = "read";
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static void usg(void) {
  fprintf(stderr, "%s mesh > off\n", me);
  exit(1);
}

int main(int argc, const char** argv) {
  USED(argc);
  int nt;
  int nv;
  int* tri;
  double* ver;

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
  off_write(nt, tri, nv, ver, stdout);
  inside_mesh_fin(tri, ver);
}
