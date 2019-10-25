#include <stdio.h>
#include <csv.h>

int
main(int argc, char **argv)
{
  int i, nr;
  double *cl;
  struct CSV *csv;

  csv = csv_read(stdin);
  nr = csv_nr(csv);
  cl = csv_field(csv, "cl");

  for (i = 0; i < nr; i++)
    printf("%g\n", cl[i]);

  csv_fin(csv);
}
