// Created by Sergey Litvinov on 31.01.2021
// Copyright 2021 ETH Zurich

void* memory_malloc(int);
void* memory_realloc(void*, int);
void memory_free(void*);
void* memory_memcpy(void*, const void*, int);

#define MALLOC(n, p)                                    \
  do {                                                  \
    *(p) = memory_malloc(n * sizeof(**(p)));            \
    if (*(p) == NULL) ERR(("alloc failed, n = %d", n)); \
  } while (0)
#define FREE(p) memory_free(p)
#define REALLOC(n, p)                                     \
  do {                                                    \
    *(p) = memory_realloc(*(p), n * sizeof(**(p)));       \
    if (*(p) == NULL) ERR(("realloc failed, n = %d", n)); \
  } while (0)
