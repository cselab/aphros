// Created by Sergey Litvinov on 24.12.2019
// Copyright 2019 ETH Zurich

#include <stdio.h>
#include <stdlib.h>
#include <aphros/march/march.h>

int main(int argc, char** argv) {
  double u[8];
  double tri[3 * 3 * MARCH_NTRI], x, y, z;
  int n, i, j, a, b, c;

  if (argv[1] && argv[1][0] == '-' && argv[1][1] == 'h') {
    fprintf(stderr, "./main -1 -1 1 1 1 1 1 -1 > OBJ\n");
    exit(2);
  }
  if (argc != 9) {
    fprintf(stderr, "needs 8 arguments (%d given)\n", argc - 1);
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
