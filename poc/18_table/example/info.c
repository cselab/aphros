#include <stdio.h>
#include <stdlib.h>
#include <table.h>

static char me[] = "table/info";

static int
cmp(const void *p, const void *q)
{
  int x, y;

  x = *(int *) p;
  y = *(int *) q;
  if (x < y)
    return -1;
  else if (x > y)
    return 1;
  else
    return 0;
}

static unsigned
hash(const void *p)
{
  int x;

  x = *(int *) p;
  return x;
}


int
main(void)
{
  struct Table *t;
  int key, val;
  int *p;

  t = table_new(0, cmp, hash);
  key = 1;
  val = 10;
  table_put(t, &key, &val);

  key = 2;
  val = 20;
  table_put(t, &key, &val);

  key = 3;
  val = 30;
  table_put(t, &key, &val);

  key = 1;
  p = table_get(t, &key);
  if (p)
    printf("%d\n", *p);
  table_free(t);
}
