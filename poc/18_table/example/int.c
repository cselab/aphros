#include <stdio.h>
#include <stdlib.h>
#include <table.h>

static char me[] = "table/info";

int
main(void)
{
  struct Table *t;
  int key, val;
  int p;

  t = table_new(0, NULL, NULL);
  key = 1;
  val = 10;
  table_put(t, &key, val);

  key = 2;
  val = 20;
  table_put(t, &key, val);

  key = 3;
  val = 30;
  table_put(t, &key, val);

  key = 20;
  table_remove(t, &key);

  key = 2;
  p = table_get(t, &key);
  if (p != NONE)
    printf("%d\n", p);
  table_free(t);
}
