enum { VTK_MAX_NF = 999 };
struct VTK {
  char *name[VTK_MAX_NF];
  double *data[VTK_MAX_NF];
  int nf, nr;
};

struct VTK *vtk_read(FILE *);
int vtk_fin(struct VTK *);
int vtk_write(struct VTK *, FILE *);
