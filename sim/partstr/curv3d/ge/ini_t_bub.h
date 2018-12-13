#define MAXBUB 10000

// file format:
// one line per bubble:
// x y z r

typedef struct {
  double x, y, z, r, ry, rz;
} Bub;

// zeroz: if true, set z=0 for all bubbles
static void read(char* fn, Bub* l, int* n, int zeroz) {
  FILE* b = fopen(fn, "r");
  if (b == NULL) {
    fprintf(stderr, "ini_t_bub.h: can't open %s\n", fn);
    abort();
  }
  int i = 0;
  double x, y, z, r, ry, rz;
  while (1) {
    int cnt = fscanf(b, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &r, &ry, &rz);
    if (cnt < 4) {
      break;
    }
    l[i].x = x;
    l[i].y = y;
    l[i].z = (zeroz ? 0. : z);
    l[i].r = r;
    l[i].ry = (cnt > 4 ? ry : r);
    l[i].rz = (cnt > 5 ? rz : r);
    ++i;
  }
  *n = i;
  g_assert(*n <= MAXBUB);

  fclose(b);

  fprintf(stderr, "Read %d bubbles\n", *n);
}

// distance to bubble interface
static double dist(Bub* b, double x, double y, double z) {
  double dx = b->x - x;
  double dy = b->y - y;
  double dz = b->z - z;
  dx /= b->r;
  dy /= b->ry;
  dz /= b->rz;
  double d = sqrt(dx * dx + dy * dy + dz * dz);
  d -= 1.;
  return d;
}

// smallest distance to bubble surface
// (negative if inside)
static double mindist(Bub* l, int n, double x, double y, double z) {
  double m = dist(&l[0], x, y, z);
  int i;
  for (i = 1; i < n; ++i) {
    double d = dist(&l[i], x, y, z);
    if (d < m) {
      m = d;
    }
  }
  return m;
}


Bub l[MAXBUB];
int n = -1;

static double ini_t(double x, double y, double z, char* fn) {
  if (n == -1) {
    read(fn, l, &n, 0);
  }
  return -mindist(l, n, x, y, z);
}

static double ini_t2(double x, double y, char* fn) {
  if (n == -1) {
    read(fn, l, &n, 1);
  }
  return -mindist(l, n, x, y, 0);
}

