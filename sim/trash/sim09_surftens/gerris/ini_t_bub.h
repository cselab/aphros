#define MAXBUB 10000
#define FBUB "../b.dat"

// file format:
// one line per bubble:
// x y z r

typedef struct {
  double x, y, z, r;
} Bub;

static void read(char* fn, Bub* l, int* n) {
  FILE* b = fopen(fn, "r");
  if (b == NULL) {
    fprintf(stderr, "ini_t_bub.h: can't open %s\n", fn);
    abort();
  }
  int i = 0;
  double x, y, z, r;
  while (fscanf(b, "%lf %lf %lf %lf\n", &x, &y, &z, &r) > 0) {
    l[i].x = x;
    l[i].y = y;
    l[i].z = z;
    l[i].r = r;
    ++i;
  }
  *n = i;
  g_assert(*n <= MAXBUB);

  fclose(b);

  fprintf(stderr, "Read %d bubbles\n", *n);
}

static double dist(Bub* b, double x, double y, double z) {
  double dx = b->x - x;
  double dy = b->y - y;
  double dz = b->z - z;
  double d = sqrt(dx * dx + dy * dy + dz * dz);
  d -= b->r;
  return d;
}

// smallest distance to bubble surface
// (negative if inside)
static double mindist(Bub* l, int n, double x, double y, double z) {
  double m = dist(&l[0], x, y, z);
  int i;
  for (i = 0; i < n; ++i) {
    double d = dist(&l[i], x, y, z);
    if (d < m) {
      m = d;
    }
  }
  return m;
}

Bub l[MAXBUB];
int n = -1;

static double ini_t(double x, double y, double z) {
  if (n == -1) {
    read(FBUB, l, &n);
  }
  return -mindist(l, n, x, y, z);
}
