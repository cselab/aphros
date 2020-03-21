#include <stdlib.h>
#include <string.h>

void* memory_malloc(int size) {
  return malloc(size);
}

void* memory_realloc(void* ptr, int size) {
  return realloc(ptr, size);
}

void memory_free(void* p) {
  free(p);
}

char* memory_strndup(const char* s, size_t n) {
  char* p = memchr(s, '\0', n);

  if (p != NULL) n = p - s;
  p = malloc(n + 1);
  if (p != NULL) {
    memcpy(p, s, n);
    p[n] = '\0';
  }
  return p;
}

void* memory_memcpy(void* a, const void* b, int n) {
  return memcpy(a, b, n);
}
