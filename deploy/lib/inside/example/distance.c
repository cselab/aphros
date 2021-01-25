#include <stdio.h>
#include <stdlib.h>
#include <inside.h>

static const char* me = "distance";
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static void usg(void) {
  fprintf(stderr, "%s [-i] -p float float float < off\n", me);
  exit(1);
}

int main(int argc, const char** argv) {
  USED(argc);
  enum { X, Y, Z };
  double p[3];
  double* ver;
  int i;
  int Invert;
  int nt;
  int nv;
  int Pflag;
  int tmp;
  int* tri;
  struct Inside* inside;

  Pflag = 0;
  Invert = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'i':
        Invert = 1;
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
  if (Invert)
    for (i = 0; i < nt; i++) {
      tmp = tri[3 * i];
      tri[3 * i] = tri[3 * i + 1];
      tri[3 * i + 1] = tmp;
    }
  inside_ini(nt, tri, ver, &inside);
  printf("%.16g\n", inside_distance(inside, p));
  inside_fin(inside);
  inside_mesh_fin(tri, ver);
}
