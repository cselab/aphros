enum { VTK_INT, VTK_FLOAT, VTK_DOUBLE };
enum { VTK_SCALAR = 1, VTK_VECTOR = 3 };
enum { VTK_CELL, VTK_POINT };

enum { VTK_MAX_NF = 99 };
struct VTK {
  double *x, *y, *z;
  int *t0, *t1, *t2;
  char* name[VTK_MAX_NF];
  int location[VTK_MAX_NF];
  int type[VTK_MAX_NF];
  int rank[VTK_MAX_NF];
  void* data[VTK_MAX_NF];
  int nv, nt, nf;
};

struct VTK* vtk_ini(
    int nv, const double* x, const double* y, const double* z, int nt,
    const int*, const int*, const int*);
struct VTK* vtk_read(FILE*);
int vtk_fin(struct VTK*);
int vtk_nv(struct VTK*);
int vtk_nt(struct VTK*);
int vtk_nf(struct VTK*);

int vtk_index(struct VTK*, const char*);
void* vtk_data(struct VTK*, const char*);

int vtk_write(struct VTK*, FILE*);
int vtk_off(struct VTK*, FILE*);
int vtk_off_color(struct VTK*, const char*, FILE*);
int vtk_remove_tri(struct VTK*, const int* flag);
int vtk_remove_orphan(struct VTK*);
int vtk_add(struct VTK*, const char*, int localtion, int type);
int vtk_remove(struct VTK*, const char*);
