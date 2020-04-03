#include "quadtree.h"

enum {NE, NW, SW, SE};
struct Node {
  double x;
  double y;
  int data;
  struct Node *elm[4];
};

struct Region {
  double L;
  double R;
  double B;
  double T;
};

static int compare(struct Node *, struct Node *);
static int found(struct Node *);
static int in_region(double, double, struct Region *region);
static int rectangle_overlaps_region(double, double, double, double, struct Region *);

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
regionsearch(struct Node *P, double L, double R, double B, double T, void found(struct Node*))
{
  double XC;
  double YC;
  struct Region region;

  region.L = L;
  region.R = R;
  region.B = B;
  region.T = T;
  XC = P->x;
  YC = P->y;
  if (in_region(XC, YC, &region)) found(P);
  if ((P->elm[NE] != NULL) && (rectangle_overlaps_region
    (XC, R, YC, T, &region)))
    regionsearch(P->elm[NE], XC, R, YC, T, found);
  if ((P->elm[NW] != NULL) && (rectangle_overlaps_region
    (L, XC, YC, T, &region)))
    regionsearch(P->elm[NW], L, XC, YC, T, found);
  if ((P->elm[SW] != NULL) && (rectangle_overlaps_region
    (L, XC, B, YC, &region)))
    regionsearch(P->elm[SW], L, XC, B, YC, found);
  if ((P->elm[SE] != NULL) && (rectangle_overlaps_region
    (XC, R, B, YC, &region)))
    regionsearch(P->elm[SE], XC, R, B, YC, found);
  return 0;
}

static int
compare(struct Node *a, struct Node *b)
{
  if (b->x > a->x)
    return b->y > a->y ? NE : SE;
  else
    return b->y > a->y ? NW : SW;
}

static int
in_region(double x, double y, struct Region *region)
{
  double LP;
  double RP;
  double TP;
  double BP;
  LP = region->L;
  RP = region->R;
  BP = region->B;
  TP = region->T;
  return LP <= x && x <= RP && BP <= y && y <= TP;
}

static int
rectangle_overlaps_region(double L, double R, double B, double T, struct Region *region)
{
  double LP;
  double RP;
  double TP;
  double BP;
  LP = region->L;
  RP = region->R;
  BP = region->B;
  TP = region->T;
  return L <= RP && R >= LP && B <= TP && T >= BP;
}

static int
found(struct Node *P)
{
  (void*)P;
  return 0;
}
