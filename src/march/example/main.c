#include <stdio.h>

#include <march.h>

int main() {
  double cube[8] = {-1, 0, 0, 0, 0, 0, 0, 1};
  double tri[3 * 3 * MARCH_NTRI];
  int n, i, j;
  int u, v, w;
  double a, b, c;

  march_cube(cube, &n, tri);
  if (n == 0) return 0;

  printf("# File type: ASCII OBJ\n");
  for (i = j = 0; i < 3 * n; i++) {
    a = tri[j++];
    b = tri[j++];
    c = tri[j++];
    printf("v %g %g %g\n", a, b, c);
  }
  for (j = 1, i = 0; i < n; i++) {
    u = j++;
    v = j++;
    w = j++;
    printf("f %d %d %d\n", u, v, w);
  }
  fprintf(stderr, "n: %d\n", n);
}
