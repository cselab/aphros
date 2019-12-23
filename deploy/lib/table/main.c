#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "err.h"
#include "memory.h"
#include "table.h"

static char me[] = "table";

struct Table {
  int size;
  int length;
  unsigned timestamp;
  struct binding {
    struct binding* link;
    int key;
    int value;
  } * *buckets;
};

static int cmp(int x, int y) {
  return x != y;
}

static unsigned hash(int x) {
  return x;
}

struct Table* table_ini(int hint) {
  struct Table* q;
  int i;

  static int primes[] = {509,  509,   1021,  2053,  4093,
                         8191, 16381, 32771, 65521, INT_MAX};
  assert(hint >= 0);
  for (i = 1; primes[i] < hint; i++)
    ;
  q = memory_malloc(sizeof(*q) + primes[i - 1] * sizeof(q->buckets[0]));
  assert(q);
  q->size = primes[i - 1];
  q->buckets = (struct binding**)(q + 1);
  for (i = 0; i < q->size; i++)
    q->buckets[i] = NULL;
  q->length = 0;
  q->timestamp = 0;
  return q;
}

int table_get(struct Table* q, int key, int* pvalue) {
  int i;
  struct binding* p;

  assert(q);
  i = hash(key) % q->size;
  for (p = q->buckets[i]; p; p = p->link)
    if (cmp(key, p->key) == 0) break;
  if (p) {
    *pvalue = p->value;
    return 0;
  } else
    return TABLE_EMPY;
}

int table_put(struct Table* q, int key, int value) {
  int i;
  struct binding* p;
  int prev;

  assert(q);
  i = hash(key) % q->size;
  for (p = q->buckets[i]; p; p = p->link)
    if (cmp(key, p->key) == 0) break;
  if (p == NULL) {
    MALLOC(1, &p);
    p->key = key;
    p->link = q->buckets[i];
    q->buckets[i] = p;
    q->length++;
    prev = TABLE_EMPY;
  } else
    prev = p->value;
  p->value = value;
  q->timestamp++;
  return prev;
}

int table_length(struct Table* q) {
  assert(q);
  return q->length;
}

int table_remove(struct Table* q, int key) {
  int i;
  struct binding** pp;

  assert(q);
  q->timestamp++;
  i = hash(key) % q->size;
  for (pp = &q->buckets[i]; *pp; pp = &(*pp)->link)
    if (cmp(key, (*pp)->key) == 0) {
      struct binding* p = *pp;
      int value = p->value;

      *pp = p->link;
      FREE(p);
      q->length--;
      return value;
    }
  return TABLE_EMPY;
}

int table_fin(struct Table* table) {
  assert(table);
  if (table->length > 0) {
    int i;
    struct binding *p, *q;

    for (i = 0; i < table->size; i++)
      for (p = table->buckets[i]; p; p = q) {
        q = p->link;
        FREE(p);
      }
  }
  FREE(table);
  return 0;
}

int* table_array(struct Table* q) {
  int i, j = 0;
  int* array;
  struct binding* p;

  assert(q);
  MALLOC(2 * q->length, &array);
  for (i = 0; i < q->size; i++)
    for (p = q->buckets[i]; p; p = p->link) {
      array[j++] = p->key;
      array[j++] = p->value;
    }
  return array;
}
