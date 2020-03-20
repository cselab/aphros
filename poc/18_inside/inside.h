struct Inside;

int inside_ini(double lo[2], double hi[2], double size, struct Inside **);
int inside_fin(struct Inside *);
int inside_update(struct Inside *, int, const int * tri, const double * x, const double * y, const double * z);
int inside_inside(struct Inside *, double, double, double);
int inside_inside_fast(struct Inside *, double, double, double);
