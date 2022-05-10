// Created by Sergey Litvinov on 31.01.2021
// Copyright 2021 ETH Zurich

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

int err_print(const char* fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  return 0;
}

void err_exit(int status) {
  exit(status);
}
