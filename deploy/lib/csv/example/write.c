#include <stdio.h>

#include <csv.h>

int main(void) {
  struct CSV* csv;

  csv = csv_read(stdin);
  csv_write_space(csv, stdout);
  csv_fin(csv);
}
