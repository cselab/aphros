#include <stdlib.h>
#include "octree.h"

enum {X, Y, Z};
enum {BNE, BNW, BSW, BSE, FNE, FNW, FSW, FSE};
struct Region {
  double lo[3];
  double hi[3];
};

static int compare(struct Node *, struct Node *);
static int in_region(double, double, double z, struct Region *region);
static int overlap(double x, double y, double z, double u, double v, double w, struct Region *);

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
  return 0;
}

int
regionsearch(struct Node *P, const double lo[3], const double hi[3], void found(struct Node*))
{
  double x;
  double y;
  double z;
  struct Region region;
  struct Node **elm;

  region.lo[X] = lo[X];
  region.lo[Y] = lo[Y];
  region.lo[Z] = lo[Z];

  region.hi[X] = hi[X];
  region.hi[Y] = hi[Y];
  region.hi[Z] = hi[Z];

  elm = P->elm;
  x = P->x;
  y = P->y;
  z = P->z;
  if (in_region(x, y, z, &region)) found(P);
#define DIR(d, x, u, y, v, z, w)				      \
  if (elm[d] != NULL && overlap(x, y, z, u, v, w, &region))	      \
    regionsearch(elm[d], lo, hi, found)				      \
  
  DIR(BNE, x, hi[X], y, hi[Y], lo[Z], z);
  DIR(BSE, x, hi[X], lo[Y], y, lo[Z], z);
  DIR(BNW, lo[X], x, y, hi[Y], lo[Z], z);
  DIR(BSW, lo[X], x, lo[Y], y, lo[Z], z);

  DIR(FNE, x, hi[X], y, hi[Y], z, hi[Z]);
  DIR(FSE, x, hi[X], lo[Y], y, z, hi[Z]);
  DIR(FNW, lo[X], x, y, hi[Y], z, hi[Z]);
  DIR(FSW, lo[X], x, lo[Y], y, z, hi[Z]);
  return 0;
}

struct Node *
node_ini(double x, double y, double z, int data)
{
  struct Node *node;
  if ((node = malloc(sizeof(*node))) == NULL)
    return NULL;
  node->x = x;
  node->y = y;
  node->z = z;
  node->data = data;
  node->elm[BNE] = node->elm[BNW] = node->elm[BSW] = node->elm[BSE] = NULL;
  node->elm[FNE] = node->elm[FNW] = node->elm[FSW] = node->elm[FSE] = NULL;
  return node;
}

static int
compare(struct Node *a, struct Node *b)
{
  if (b->z > a->z) {
    if (b->x > a->x)
      return b->y > a->y ? FNE : FSE;
    else
      return b->y > a->y ? FNW : FSW;
  } else {
    if (b->x > a->x)
      return b->y > a->y ? BNE : BSE;
    else
      return b->y > a->y ? BNW : BSW;    
  }
}

static int
in_region(double x, double y, double z, struct Region *region)
{
  const double *lo;
  const double *hi;
  lo = region->lo;
  hi = region->hi;
  return
    lo[X] <= x && x <= hi[X] &&
    lo[Y] <= y && y <= hi[Y] &&
    lo[Z] <= z && z <= hi[Z];
}

static int
overlap(double x, double y, double z, double u, double v, double w, struct Region *region)
{
  const double *lo;
  const double *hi;
  lo = region->lo;
  hi = region->hi;
  return
    x <= hi[X] && u >= lo[X] &&
    y <= hi[Y] && v >= lo[Y] &&
    z <= hi[Z] && w >= lo[Z];
}
