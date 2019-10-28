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

  t = table_ini(0, NULL, NULL);
  key = 40;
  val = 1;
  table_put(t, &key, val);

  key = 549;
  val = 2;
  table_put(t, &key, val);

  key = 40;
  p = table_get(t, &key);
  if (p != TABLE_EMPY)
    printf("%d\n", p);
  table_fin(t);
}
