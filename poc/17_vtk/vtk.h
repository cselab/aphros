enum { VTK_MAX_NF = 999 };
enum { VTK_INT, VTK_FLOAT };
enum { VTK_SCALAR = 1, VTK_VECTOR = 3 };

struct VTK {
  double *x, *y, *z;
  int *t0, *t1, *t2;
  char *name[VTK_MAX_NF];
  int type[VTK_MAX_NF];
  int rank[VTK_MAX_NF];
  void *data[VTK_MAX_NF];
  int nv, nt, nf;
};

struct VTK *vtk_read(FILE *);
int vtk_fin(struct VTK *);
int vtk_nv(struct VTK *);
int vtk_nt(struct VTK *);
int vtk_write(struct VTK *, FILE *);
