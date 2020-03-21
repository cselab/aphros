#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <inside.h>

const char* me = "off";

static void usg(void) {
  fprintf(stderr, "%s < off > off\n", me);
  exit(1);
}

int
main(int argc, const char** argv) {
  int nt;
  int nv;
  int *tri;
  double *ver;
  int status;

  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    default:
      fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
      exit(2);
    }
  off_read(stdin, &status, &nt, &tri, &nv, &ver);
  off_write(nt, tri, nv, ver, stdout);
  off_fin(tri, ver);
}
