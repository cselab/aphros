typedef struct {
  double x, y, z;
} coord;
typedef struct { int i; } scalar;
typedef struct {
  int i;
  int j;
  int k;
  int level;
} Point;

void init_grid(int n);
double line_alpha(double, coord);
double embed_interpolate (Point point, scalar s, coord p);
