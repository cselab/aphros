enum { VTK_MAX_NF = 999 };
struct VTK {
  char *name[VTK_MAX_NF];
  int *type[VTK_MAX_NF];
  int *rank[VTK_MAX_NF];
  void *data[VTK_MAX_NF];
  int nf, nr;
};

struct VTK *vtk_read(FILE *);
int vtk_fin(struct VTK *);
int vtk_write(struct VTK *, FILE *);
