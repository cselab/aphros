#include <stdio.h>

#include <csv.h>

int main(void) {
  struct CSV* csv;
  double *age, *x;
  int i, n;

  csv = csv_ini(10);

  csv_add(csv, "age");
  csv_add(csv, "x");
  age = csv_field(csv, "age");
  x = csv_field(csv, "x");
  n = csv_nr(csv);
  for (i = 0; i < n; i++) {
    age[i] = i;
    x[i] = 10 * i;
  }

  csv_write(csv, stdout);
  csv_fin(csv);
}
