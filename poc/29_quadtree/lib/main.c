#include <stdio.h>
#include "quadtree.h"

enum {NE, SE, NW, SW};
struct Node {
  double x;
  double y;
  int data;
  struct Node *elm[4];
};
static int compare(struct Node *, struct Node *);

int
insert(struct Node *K, struct Node *R)
{
  int direction;
  direction = compare(R, K);
  while (R->elm[direction] != NULL) {
    R = R->elm[direction];
    direction = compare(R, K);
  }
  R->elm[direction] = K;
}

static int
compare(struct Node *a, struct Node *b)
{
  if (b->x > a->x)
    return b->y > a->y ? NE : SE;
  else
    return b->y > a->y ? NW : SW;
}
