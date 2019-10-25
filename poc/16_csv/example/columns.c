#include <stdio.h>
#include <stdlib.h>
#include <csv.h>

static const char me[] = "column";

#define FMT "%.20g"

int
main(int argc, char **argv)
{
  int i, j, nf, nr;
  char *name;
  double *field[99], *f;
  struct CSV *csv;

  argv++;

  csv = csv_read(stdin);
  nr = csv_nr(csv);
  nf = 0;
  for (;;) {
    name = argv[0];
    argv++;
    if (name == NULL)
      break;
    f = csv_field(csv, name);
    if (f == NULL) {
      fprintf(stderr, "%s: unknow filed '%s'\n", me, name);
      exit(2);
    }
    field[nf++] = f;
  }

  for (i = 0; i < nr; i++) {
    for (j = 0; j < nf; j++) {
      if (j > 0)
        fputc(' ', stdout);
      fprintf(stdout, FMT, field[j][i]);
    }
    fputc('\n', stdout);
  }


  csv_fin(csv);
}
