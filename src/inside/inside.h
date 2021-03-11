// Created by Sergey Litvinov on 31.01.2021
// Copyright 2021 ETH Zurich

#ifdef __cplusplus
extern "C" {
#endif
struct Inside;
struct InsideInfo {
  double size;
  int nx;
  int ny;
  int max_tri;
  int min_tri;
};
int inside_ini(int, const int* tri, const double* ver, struct Inside**);
int inside_fin(struct Inside*);
int inside_inside(struct Inside*, const double r[3]);
int inside_inside_naive(struct Inside*, const double r[3]);
double inside_distance(struct Inside*, const double[3]);
double inside_distance_naive(struct Inside*, const double[3]);
int inside_box(struct Inside*, double[3], double[3]);
int inside_info(struct Inside*, struct InsideInfo*);
int inside_fwrite(struct Inside*, FILE*);
int inside_mesh_read(const char*, int* nt, int** tri, int* nv, double** ver);
int inside_mesh_fin(int* tri, double* ver);
int off_write(int nt, const int*, int nv, const double*, FILE*);
int ply_write(int nt, const int*, int nv, const double*, FILE*);
int stl_write(int nt, const int*, int nv, const double*, FILE*);
#ifdef __cplusplus
}
#endif
