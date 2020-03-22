struct Inside;

int inside_ini(int, const int * tri, const double * ver, struct Inside **);
int inside_fin(struct Inside *);
int inside_inside(struct Inside *, const double r[3]);

int inside_mesh_read(FILE *, int * nt, int ** tri, int * nv, double ** ver);
int inside_mesh_fin(int *tri, double *ver);

int off_write(int nt, const int *, int nv, const double *, FILE *);
