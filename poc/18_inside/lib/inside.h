struct Inside;

int inside_ini(struct Inside **);
int inside_fin(struct Inside *);
int inside_update(struct Inside *, int, const int * tri, const double * x, const double * y, const double * z);
int inside_inside(struct Inside *, double, double, double);

int off_read(FILE *, int * status, int * nt, int ** tri, int * nv, double ** x, double ** y, double ** z);
int off_fin(int *tri, double *, double *, double *);
int off_write(int nt, int *, int nv, double *, double *, double *, FILE *);

