#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <table.h>

static char me[] = "table/stream";
enum { N = 999 };
static int eq(const char*, const char*);
static char* get(char*, FILE*);

#define ARG(x)                                            \
  do {                                                    \
    c = strtok(NULL, sep);                                \
    if (c == NULL) {                                      \
      fprintf(stderr, "%s: wrong command '%s'\n", me, p); \
      exit(2);                                            \
    }                                                     \
    (x) = atoi(c);                                        \
  } while (0)

int main(void) {
  char s[N], p[N], sep[] = " \t\n";
  char* c;
  static struct Table* t;
  int i, j, x, y, n, *a;

  t = table_ini(0);
  while (get(s, stdin)) {
    strncpy(p, s, N);
    if (c = strtok(s, sep)) {
      if (eq(c, "put")) {
        ARG(x);
        ARG(y);
        table_put(t, x, y);
      } else if (eq(c, "get")) {
        ARG(x);
        if (table_get(t, x, &y) == TABLE_EMPY)
          printf("empty\n");
        else
          printf("%d\n", y);
      } else if (eq(c, "remove")) {
        ARG(x);
        table_remove(t, x);
      } else if (eq(c, "length")) {
        printf("%d\n", table_length(t));
      } else if (eq(c, "array")) {
        a = table_array(t);
        n = table_length(t);
        for (i = j = 0; i < n; i++) {
          printf("[%d %d]", a[j++], a[j++]);
        }
        if (n > 0) putc('\n', stdout);
        free(a);
      } else {
        fprintf(stderr, "%s: unknown command '%s'\n", me, c);
        exit(2);
      }
    }
  }
  table_fin(t);
}

static int eq(const char* a, const char* b) {
  return strncmp(a, b, N) == 0;
}

static char* get(char* s, FILE* f) {
  int n;

  if (fgets(s, N, f) == NULL) return NULL;
  n = strlen(s);
  if (s[n - 1] == '\n') s[n - 1] = '\0';
  return s;
}
