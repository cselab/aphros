struct Bbox {
  double lo[3];
  double hi[3];
};

int bbox_ini(struct Bbox** pq) {
  struct Bbox* q;

  MALLOC(1, &q);
  *pq = q;
  return 0;
}

int bbox_fin(struct Bbox* q) {
  FREE(q);
  return 0;
}

int bbox_update(struct Bbox* q, int n, const double* ver) {
  int i;
  int j;
  double x;
  double y;
  double z;
  for (j = i = 0; i < n; i++) {
    x = ver[j++];
    y = ver[j++];
    z = ver[j++];
    if (i == 0) {
      q->lo[X] = q->hi[X] = x;
      q->lo[Y] = q->hi[Y] = y;
      q->lo[Z] = q->hi[Z] = z;
    } else {
      if (x < q->lo[X]) q->lo[X] = x;
      if (y < q->lo[Y]) q->lo[Y] = y;
      if (z < q->lo[Z]) q->lo[Z] = z;
      if (x > q->hi[X]) q->hi[X] = x;
      if (y > q->hi[Y]) q->hi[Y] = y;
      if (z > q->hi[Z]) q->hi[Z] = z;
    }
  }
  return 0;
}

int bbox_inside(struct Bbox* q, const double r[3]) {
#define CM(D) (lo[D] < r[D] && r[D] < hi[D])
  double *lo, *hi;

  lo = q->lo;
  hi = q->hi;
  return CM(X) && CM(Y) && CM(Z);
}

double bbox_zhi(struct Bbox* q) {
  return q->hi[Z];
}

int bbox_lo(struct Bbox* q, const double** x) {
  *x = q->lo;
  return 0;
}

int bbox_hi(struct Bbox* q, const double** x) {
  *x = q->hi;
  return 0;
}
