#include <stdio.h>
#include <stdlib.h>

#include <march.h>

int main(int argc, char** argv) {
  double u[8];
  double tri[3 * 3 * MARCH_NTRI], x, y, z;
  int n, i, j, a, b, c;

  if (argc != 9) {
    fprintf(stderr, "needs eith arguments (%d given)\n", argc - 1);
    exit(2);
  }
  for (i = 0; i < 8; i++)
    u[i] = atof(*++argv);
  march_cube(u, &n, tri);
  printf("# File type: ASCII OBJ\n");
  for (i = j = 0; i < 3 * n; i++) {
    x = tri[j++];
    y = tri[j++];
    z = tri[j++];
    printf("v %g %g %g\n", x, y, z);
  }
  for (j = 1, i = 0; i < n; i++) {
    a = j++;
    b = j++;
    c = j++;
    printf("f %d %d %d\n", a, b, c);
  }
}
