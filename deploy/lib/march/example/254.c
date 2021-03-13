#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <march.h>

static const char* me = "254";

static void tobin(int x, double* b) {
  int i;
  for (i = 0; i < 8; i++) {
    b[i] = x % 2 == 0 ? -1 : 1;
    x >>= 1;
  }
}

int main() {
  double cube[8];
  double tri[3 * 3 * MARCH_NTRI];
  char name[999];
  int id;
  int n, i, j;
  int u, v, w;
  double a, b, c;
  FILE* file;

  for (id = 1; id < 255; id++) {
    tobin(id, cube);
    sprintf(name, "a.%03d.obj", id);
    if ((file = fopen(name, "w")) == NULL) {
      fprintf(stderr, "%s: fail to write to '%s'\n", me, name);
      exit(2);
    }
    march_cube(cube, &n, tri);
    assert(n != 0);
    fprintf(file, "# File type: ASCII OBJ\n");
    for (i = j = 0; i < 3 * n; i++) {
      a = tri[j++];
      b = tri[j++];
      c = tri[j++];
      fprintf(file, "v %g %g %g\n", a, b, c);
    }
    for (j = 1, i = 0; i < n; i++) {
      u = j++;
      v = j++;
      w = j++;
      fprintf(file, "f %d %d %d\n", u, v, w);
    }
    fclose(file);
  }
}
