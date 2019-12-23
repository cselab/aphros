#include <stdio.h>

#include <csv.h>

int main(void) {
  struct CSV* csv;
  double* age;
  int i, n;

  csv = csv_read(stdin);

  csv_add(csv, "age");
  age = csv_field(csv, "age");
  n = csv_nr(csv);
  for (i = 0; i < n; i++)
    age[i] = i;

  csv_write_space(csv, stdout);
  csv_fin(csv);
}
