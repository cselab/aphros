#include <stdio.h>
#include <stdlib.h>

#include <csv.h>

static const char me[] = "column";

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

int main(int argc, char** argv) {
  int i, nr;
  double* field;
  struct CSV* csv;

  USED(argc);
  argv++;

  csv = csv_read(stdin);
  nr = csv_nr(csv);
  if (argv[0] == NULL) {
    fprintf(stderr, "%s: needs an argument\n", me);
    exit(2);
  }
  field = csv_field(csv, argv[0]);
  if (field != NULL)
    for (i = 0; i < nr; i++)
      printf("%.20g\n", field[i]);
  csv_fin(csv);
}
