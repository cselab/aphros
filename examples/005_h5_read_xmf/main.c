#include <stdio.h>
#include <stdlib.h>
#include <h5serial.h>

int main(int argc, char** argv) {
  enum { X, Y, Z };
  const char *path, *me;
  char name[999];
  double ori[3], spa;
  int siz[3];

  me = argv[0];
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        fprintf(stderr, "%s PREFIX\n", me);
        exit(2);
        break;
      default:
        fprintf(stderr, "%s: unknow option\n", me);
        exit(2);
    }
  if ((path = argv++ [0]) == NULL) {
    fprintf(stderr, "%s: missing xmf path\n", me);
    exit(2);
  }
  if (h5_read_xmf(path, name, ori, &spa, siz) != 0) {
    fprintf(stderr, "%s: fail to read '%s'\n", me, path);
    exit(2);
  }
  printf("name: %s\n", name);
  printf("ori: %g %g %g\n", ori[X], ori[Y], ori[Z]);
  printf("spa: %g\n", spa);
  printf("siz: %d %d %d\n", siz[X], siz[Y], siz[Z]);
}
