#include <stdio.h>
#include <stdlib.h>
#include <inside.h>

static const char* me = "point";
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
  double p[3];
  int Pflag;
  struct Inside* inside;

  Pflag = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'p':
        argv++;
        if (argv[0] == NULL || argv[1] == NULL || argv[2] == NULL) {
          fprintf(stderr, "%s: -p needs three numbers\n", me);
          exit(2);
        }
        p[X] = atof(*argv++);
        p[Y] = atof(*argv++);
        p[Z] = atof(*argv);
        Pflag = 1;
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(2);
    }
  if (Pflag == 0) {
    fprintf(stderr, "%s: needs -p option\n", me);
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
  printf("%d\n", inside_inside(inside, p));
  inside_fin(inside);
  inside_mesh_fin(tri, ver);
}
