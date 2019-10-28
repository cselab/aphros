#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include "err.h"
#include "memory.h"
#include "table.h"

static char me[] = "table";

struct Table {
  int size;
  int (*cmp)(const void *x, const void *y);
  unsigned (*hash)(const void *key);
  int length;
  unsigned timestamp;
  struct binding {
    struct binding *link;
    const void *key;
    int value;
  } **buckets;
};

static int
cmpint(const void *p, const void *q)
{
  int x, y;

  x = *(int *) p;
  y = *(int *) q;
  return x != y;
}

static unsigned
hashint(const void *p)
{
  int x;

  x = *(int *) p;
  return x;
}

struct Table *
table_ini(int hint,
          int cmp(const void *x, const void *y),
          unsigned hash(const void *key))
{
  struct Table *q;
  int i;

  static int primes[] = { 509, 509, 1021, 2053, 4093,
    8191, 16381, 32771, 65521, INT_MAX
  };
  assert(hint >= 0);
  for (i = 1; primes[i] < hint; i++);
  q = memory_malloc(sizeof(*q) + primes[i - 1] * sizeof(q->buckets[0]));
  assert(q);
  q->size = primes[i - 1];
  q->cmp = cmp ? cmp : cmpint;
  q->hash = hash ? hash : hashint;
  q->buckets = (struct binding **) (q + 1);
  for (i = 0; i < q->size; i++)
    q->buckets[i] = NULL;
  q->length = 0;
  q->timestamp = 0;
  return q;
}

int
table_get(struct Table *q, const void *key)
{
  int i;
  struct binding *p;

  assert(q);
  assert(key);
  i = (*q->hash) (key) % q->size;
  for (p = q->buckets[i]; p; p = p->link)
    if ((*q->cmp) (key, p->key) == 0)
      break;
  return p ? p->value : TABLE_EMPY;
}

int
table_put(struct Table *q, const void *key, int value)
{
  int i;
  struct binding *p;
  int prev;

  assert(q);
  assert(key);
  i = (*q->hash) (key) % q->size;
  for (p = q->buckets[i]; p; p = p->link)
    if ((*q->cmp) (key, p->key) == 0)
      break;
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

int
table_length(struct Table *q)
{
  assert(q);
  return q->length;
}

void
table_map(struct Table *q,
          void apply(const void *key, int *value, void *cl), void *cl)
{
  int i;
  unsigned stamp;
  struct binding *p;

  assert(q);
  assert(apply);
  stamp = q->timestamp;
  for (i = 0; i < q->size; i++)
    for (p = q->buckets[i]; p; p = p->link) {
      apply(p->key, &p->value, cl);
      assert(q->timestamp == stamp);
    }
}

int
table_remove(struct Table *q, const void *key)
{
  int i;
  struct binding **pp;

  assert(q);
  assert(key);
  q->timestamp++;
  i = (*q->hash) (key) % q->size;
  for (pp = &q->buckets[i]; *pp; pp = &(*pp)->link)
    if ((*q->cmp) (key, (*pp)->key) == 0) {
      struct binding *p = *pp;
      int value = p->value;

      *pp = p->link;
      FREE(p);
      q->length--;
      return value;
    }
  return TABLE_EMPY;
}

void
table_fin(struct Table *table)
{
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
}
