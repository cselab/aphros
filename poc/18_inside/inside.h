struct Inside;
struct He;

int surface_ini(double lo[2], double hi[2], double size, struct Inside **);
int surface_fin(struct Inside *);
int surface_update(struct Inside *, struct He *, const double * x, const double * y, const double * z);
int surface_inside(struct Inside *, double, double, double);
int surface_inside_fast(struct Inside *, double, double, double);
