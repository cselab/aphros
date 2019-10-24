#include <stdio.h>
#include "h.h"

#define N (1024)

int
main()
{
  int k;

  h_ini(10);
  h_enter(1.0, 42);
  h_enter(45.0, 56);

  k = h_find(1.0);
  fprintf(stderr, "k: %d\n", k == -1);
  h_fin();
}
