#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "err.h"

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
