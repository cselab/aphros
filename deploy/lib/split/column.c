#include <stdio.h>
#include <stdlib.h>

#include <csv.h>

static const char me[] = "column";

#define FMT "%.20g"
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

static void usg() {
  fprintf(stderr, "%s [-s separator] [-k] [field..] < CSV\n", me);
  exit(1);
}

int main(int argc, char** argv) {
  int i, j, nf, nr;
  int Key;
  char* name;
  double* field[99];
  double* f;
  char* key[99];
  struct CSV* csv;
  char* Sep = " ";

  USED(argc);
  Key = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 's':
        argv++;
        Sep = argv[0];
        break;
      case 'k':
        Key = 1;
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(1);
    }
  if (argv[0] == NULL) {
    fprintf(stderr, "%s: needs an argument\n", me);
    exit(2);
  }
  csv = csv_read(stdin);
  nr = csv_nr(csv);
  nf = 0;
  for (;;) {
    name = argv[0];
    argv++;
    if (name == NULL) break;
    f = csv_field(csv, name);
    if (f == NULL) {
      fprintf(stderr, "%s: unknown filed '%s'\n", me, name);
      exit(2);
    }
    key[nf] = name;
    field[nf] = f;
    nf++;
  }

  if (Key) {
    for (j = 0; j < nf; j++) {
      if (j > 0) fputs(Sep, stdout);
      fputs(key[j], stdout);
    }
    fputc('\n', stdout);
  }

  for (i = 0; i < nr; i++) {
    for (j = 0; j < nf; j++) {
      if (j > 0) fputs(Sep, stdout);
      fprintf(stdout, FMT, field[j][i]);
    }
    fputc('\n', stdout);
  }
  csv_fin(csv);
}
