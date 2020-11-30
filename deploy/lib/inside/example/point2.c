#include <stdio.h>
#include <stdlib.h>
#include <inside.h>

static const char* me = "point2";
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static void usg(void) {
  fprintf(stderr, "%s -p float float float < off\n", me);
  exit(1);
}

int main(int argc, const char** argv) {
  USED(argc);
  enum { X, Y, Z };
  int nt;
  int nv;
  int* tri;
  double* ver;
  int data;
  int in;
  int out;
  double p[3];
  struct Inside* inside;
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
    fprintf(stderr, "%s: fail to read mesh: '%s'\n", me, argv[0]);
    exit(2);
  }
  inside_ini(nt, tri, ver, &inside);

  in = out = 0;
  while (scanf("%lf %lf %lf", &p[X], &p[Y], &p[Z]) == 3) {
    data = inside_inside(inside, p);
    if (data == 0)
      out++;
    else
      in++;
    printf("%d\n", data);
  }
  if (getenv("LOG") != NULL) {
    fprintf(stderr, "in:  %08d\n", in);
    fprintf(stderr, "out: %08d\n", out);
  }
  inside_fin(inside);
  inside_mesh_fin(tri, ver);
}
