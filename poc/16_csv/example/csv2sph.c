#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <csv.h>
#include <table.h>

enum { N = 999 };
static const char me[] = "csv2sph";
#include "util.h"

#define GET(f, r)							\
  if ((*r = csv_field(csv, f)) == NULL) {				\
    fprintf(stderr, "%s: no field '%s' in '%s'\n",			\
	    me, (f), *argv);						\
    exit(2);								\
  }


#define	USED(x)		if(x);else{}
#define MALLOC(n, p)							\
    do {								\
	*(p) = malloc((n)*sizeof(**(p)));				\
	if (*(p) == NULL)  {						\
	    fprintf(stderr, "%s: alloc failed, n = %d", me, n);		\
	    exit(2);							\
	}								\
    } while(0)
#define REALLOC(n, p)							\
    do {								\
      *(p) = realloc(*(p), (n)*sizeof(**(p)));				\
      if (*(p) == NULL)  {						\
	fprintf(stderr, "%s: realloc failed, n = %d", me, n);		\
	    exit(2);							\
	}								\
    } while(0)

static void
usg()
{
  fprintf(stderr, "%s -p prefix [csv..]\n", me);
  exit(1);
}

int
main(int argc, char **argv)
{
  char *Prefix;
  char output[N];
  FILE *file;
  struct CSV *csv;
  double *x;
  double *y;
  double *z;
  double *r;

  USED(argc);
  Prefix = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    case 'p':
      argv++;
      if ((Prefix = *argv) == NULL) {
        fprintf(stderr, "%s: -p needs an argument\n", me);
        exit(2);
      }
      break;
    default:
      fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
      exit(1);
    }
  if (Prefix == NULL) {
    fprintf(stderr, "%s: prefix (-p) is not given\n", me);
    exit(1);
  }
  if (*argv == NULL) {
    fprintf(stderr, "%s: csv file is not given\n", me);
    exit(1);
  }

  for ( ; *argv != NULL; argv++) {
    if ((file = fopen(*argv, "r")) == NULL) {
      fprintf(stderr, "%s: fail to open '%s'\n", me, *argv);
      exit(1);
    }
    if ((csv = csv_read(file)) == NULL) {
      fprintf(stderr, "%s: not a cvs file '%s'\n", me, *argv);
      exit(1);
    }
    fclose(file);

    GET("x", &x);
    GET("y", &y);
    GET("z", &z);
    GET("r", &r);
    if (util_name(Prefix, *argv, output) != 0) {
      fprintf(stderr, "%s: util_name failed\n", me);
      exit(2);
    }
    if ((file = fopen(output, "w")) == NULL) {
      fprintf(stderr, "%s: fail to write to '%s'\n", me, output);
      exit(2);
    }
    fclose(file);
  }
}
