#include <stdio.h>
#include <stdlib.h>
#include <inside.h>

const char* me = "point2";

static void usg(void) {
  fprintf(stderr, "%s -p float float float < off\n", me);
  exit(1);
}

int
main(int argc, const char** argv) {
  enum {X, Y, Z};
  int nt;
  int nv;
  int *tri;
  double *ver;
  int status;
  int data;
  int in;
  int out;
  double p[3];
  FILE *file;
  struct Inside *inside;
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
    fprintf(stderr, "%s: off file is missing\n", me);
    exit(2);
  }
  if ((file = fopen(argv[0], "r")) == NULL) {
    fprintf(stderr, "%s: failed to open '%s'\n", me, argv[0]);
    exit(2);    
  }
  if (off_read(file, &status, &nt, &tri, &nv, &ver) != 0) {
    fprintf(stderr, "%s: off_read failed\n", me);
    exit(2);
  }
  fclose(file);
  if (status != 0) {
    fprintf(stderr, "%s: not an off file\n", me);
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
    fprintf(stderr, "in:  %8.d\n", in);
    fprintf(stderr, "out: %8.d\n", out);
  }
  inside_fin(inside);
  off_fin(tri, ver);
}
