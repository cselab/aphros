#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "err.h"
#include "line.h"
#include "memory.h"

enum { N = 999 };
enum { LINE_N = 999 };
static struct {
  char* buf[LINE_N];
  int i;
  int n;
} line_info;
static int line(char*, FILE*);
static int line0(char*, FILE*);

int line_ini(void) {
  line_info.n = 0;
  line_info.i = 0;
  return 0;
}

int line_get(char* p, FILE* f) {
  int n, i, status;
  char** buf;
  char s[N];

  n = line_info.n;
  i = line_info.i;
  buf = line_info.buf;

  if (i >= n) {
    status = line(s, f);
    if (status != 0) goto fail;
    if (n >= LINE_N) {
      MSG(("n=%d >= LINE_N", n));
      goto fail;
    }
    buf[n] = memory_strndup(s, N);
    if (buf[n] == NULL) {
      MSG(("memory_strndup failed, n = %d", n));
      goto fail;
    }
    strncpy(p, buf[n++], N);
    i++;
  } else
    strncpy(p, buf[i++], N);
  line_info.n = n;
  line_info.i = i;
  return 0;
fail:
  return 1;
}

int line_unget(void) {
  line_info.i--;
  if (line_info.i < 0) return 1;
  return 0;
}

int line_fin(void) {
  int n, i;
  char** buf;

  n = line_info.n;
  buf = line_info.buf;
  for (i = 0; i < n; i++)
    FREE(buf[i]);
  return 0;
}

int line_write(FILE* f) {
  int n, i;
  char** buf;

  n = line_info.n;
  buf = line_info.buf;

  for (i = 0; i < n; i++) {
    fputs(buf[i], f);
    putc('\n', f);
  }
  return 0;
}

static int line(char* s, FILE* f) {
  do
    if (line0(s, f) != 0) return 1;
  while (s[0] == '\0');
  return 0;
}

static int line0(char* s, FILE* f) {
  int n;

  if (fgets(s, N, f) == NULL) return 1;
  n = strlen(s);
  if (n > 0 && s[n - 1] == '\n') s[n - 1] = '\0';
  return 0;
}
