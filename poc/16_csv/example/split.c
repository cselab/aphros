#include <stdio.h>
#include <stdlib.h>
#include <csv.h>
#include <table.h>

static const char me[] = "split";

#define	USED(x)		if(x);else{}

static void
usg()
{
  fprintf(stderr, "%s -f field [csv..]\n", me);
  exit(1);
}

int
main(int argc, char **argv)
{
  int i, nr;
  double *field;
  const char *name;
  struct CSV *csv;

  USED(argc);
  name = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    case 'f':
      argv++;
      name = argv[0];
      break;
    default:
      fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
      exit(1);
    }
  if (name == NULL) {
    fprintf(stderr, "%s: -f is not set\n", me);
    exit(1);
  }
  if (*argv == NULL) {
    fprintf(stderr, "%s: csv file is not given\n", me);
    exit(1);
  }

  FILE *file;
  struct Table *t;

  if ((file = fopen(*argv, "r")) == NULL) {
    fprintf(stderr, "%s: fail to open '%s'\n", me, *argv);
    exit(1);
  }
  if ((csv = csv_read(file)) == NULL) {
    fprintf(stderr, "%s: not a cvs file '%s'\n", me, *argv);
    exit(1);
  }
  fclose(file);
  if ((field = csv_field(csv, name)) == NULL) {
    fprintf(stderr, "%s: no field '%s' in file '%s'\n", me, name, *argv);
    exit(1);
  }
  nr = csv_nr(csv);
  t = table_ini(9999);
  for (i = 0; i < nr; i++) {
    table_put(t, field[i], i);
  }
  table_fin(t);
  csv_fin(csv);
}
