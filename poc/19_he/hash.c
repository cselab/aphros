#include <stdio.h>
#include "co/memory.h"
#include "co/err.h"
#include "co/hash.h"

#define T HeHash

typedef struct Node Node;
struct Node {
    int j;
    int v;                      /* value */
    Node *nxt;
};
struct T {
    int n;
    Node **node;
};

int
he_hash_ini(int n, T ** pq)
{
    T *q;
    int i;

    MALLOC(1, &q);
    MALLOC(q->n = n, &q->node);
    for (i = 0; i < n; i++)
        q->node[i] = NULL;
    *pq = q;
    return CO_OK;
}

static void
free_node(Node * p)
{
    Node *n;

    while (p != NULL) {
        n = p->nxt;
        FREE(p);
        p = n;
    }
}

int
he_hash_fin(T * q)
{
    int n, i;

    n = q->n;
    for (i = 0; i < n; i++)
        free_node(q->node[i]);
    FREE(q->node);
    FREE(q);
    return CO_OK;
}

int
he_hash_set(T * q, int i, int j, int v)
{
    int n;
    Node *nxt, *prv, *new;

    n = q->n;
    if (0 > i || i >= n || 0 > j || j >= n)
        ERR(CO_INDEX, "[%d %d], n = %d", i, j, n);

    prv = q->node[i];
    while (prv != NULL) {
        if (prv->j == j) {
            prv->v = v;
            return CO_OK;
        }
        nxt = prv->nxt;
        if (nxt == NULL)
            break;
        prv = nxt;
    }

    MALLOC(1, &new);
    new->j = j;
    new->v = v;
    new->nxt = NULL;
    if (prv == NULL)
        q->node[i] = new;
    else
        prv->nxt = new;
    return CO_OK;
}

int
he_hash_get(T * q, int i, int j)
{
    int n;
    Node *node;

    n = q->n;
    if (0 > i || i >= n || 0 > j || j >= n)
        ERR(CO_INDEX, "[%d %d], n = %d", i, j, n);
    node = q->node[i];
    while (node != NULL) {
        if (node->j == j)
            return node->v;
        node = node->nxt;
    }
    return -1;
}
