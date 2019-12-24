#include <stdio.h>
#include <stdlib.h>
#include <h5serial.h>

int main(int argc, char** argv) {
  enum { X, Y, Z };
  const char *path;
  char name[999];
  double ori[3], spa;
  int siz[3];
  if ((path = argv[1]) == NULL) {
    fprintf(stderr, "missing xmf path\n");
    exit(2);
  }
  if (h5_read_xmf(path, name, ori, &spa, siz) != 0) {
    fprintf(stderr, "fail to read '%s'\n", path);
    exit(2);
  }
  printf("name: %s\n", name);
  printf("ori: %g %g %g\n", ori[X], ori[Y], ori[Z]);
  printf("spa: %g\n", spa);
  printf("siz: %d %d %d\n", siz[X], siz[Y], siz[Z]);
}
