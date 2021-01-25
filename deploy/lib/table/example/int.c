#include <stdio.h>
#include <stdlib.h>

#include <table.h>

static char me[] = "table/int";

int main(void) {
  struct Table* t;
  int key, val, status;
  int p;

  t = table_ini(0);
  key = 40;
  val = 18;
  table_put(t, key, val);

  key = 549;
  val = 2;
  table_put(t, key, val);

  key = 40;
  status = table_get(t, key, &p);
  if (status != TABLE_EMPY) printf("%d\n", p);
  table_fin(t);
}
