// Created by Sergey Litvinov on 24.12.2019
// Copyright 2019 ETH Zurich

#include <stdio.h>
#include <stdlib.h>
#include <h5serial.h>

int main(int argc, char** argv) {
  enum { X, Y, Z };
  char name[999];
  double ori[3], spa;
  int siz[3];
  if (h5_read_xmf(argv[1], name, ori, &spa, siz) != 0) {
    fprintf(stderr, "fail to read '%s'\n", argv[1]);
    exit(2);
  }
  printf("name: %s\n", name);
  printf("ori: %g %g %g\n", ori[X], ori[Y], ori[Z]);
  printf("spa: %g\n", spa);
  printf("siz: %d %d %d\n", siz[X], siz[Y], siz[Z]);
}
