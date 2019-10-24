#include <stdlib.h>

void *
memory_malloc(int size)
{
  return malloc(size);
}

void
memory_free(void *p)
{
  free(p);
}

void*
memory_realloc(void *ptr, size_t size)
{
    return realloc(ptr, size);
}
