#include <stdio.h>

int foo(void f(int, int), int a, int b) {
  int i;
  for (i = 0; i < 10; i++) {
    printf("preved\n");
    f(a, b);
  }
}
