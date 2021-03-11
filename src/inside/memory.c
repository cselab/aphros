// Created by Sergey Litvinov on 31.01.2021
// Copyright 2021 ETH Zurich

#include <stdlib.h>
#include <string.h>
#include "memory.h"

void* memory_malloc(int size) {
  return malloc(size);
}

void* memory_realloc(void* ptr, int size) {
  return realloc(ptr, size);
}

void memory_free(void* p) {
  free(p);
}

void* memory_memcpy(void* a, const void* b, int n) {
  return memcpy(a, b, n);
}
