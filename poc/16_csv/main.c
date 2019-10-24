#include <stdio.h>
#include <string.h>
#include "err.h"
#include "memory.h"
#include "csv.h"

char me[] = "csv";
enum { N = 9999 };

#define LINE(s, f)				\
    if (line(s, f) != 0)			\
	ERR(("fail to read"));
static int line(char *, FILE *);

int
csv_read(FILE * f, struct CSV *q)
{
  char s[N], *name;
  int nf;
  LINE(s, f);

  nf = 0;
  name = strtok(s, ",");
  while (name != NULL) {
      nf++;
      MSG(("%s", name));
      name = strtok(NULL, ",");
  }
  q->nf = nf;
  return 0;
}

static int
line(char *s, FILE * f)
{
  int n;

  if (fgets(s, N, f) == NULL)
    return 1;
  n = strlen(s);
  if (n > 0 && s[n - 1] == '\n')
    s[n - 1] = '\0';
  return 0;
}
