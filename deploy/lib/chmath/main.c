#include <stdio.h>
#include <string.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "chmath.h"
#include "err.h"

static char* me = "chmath";

static double get(const gsl_matrix* m, int i, int j) {
  return gsl_matrix_get(m, i, j);
}

int chmath_eig_values(const double A[6], /**/ double VAL[3]) {
  enum { XX, XY, XZ, YY, YZ, ZZ };
  enum { YX = XY, ZX = XZ, ZY = YZ };
  enum { X, Y, Z };

  double B[3 * 3];
  gsl_matrix_view m;
  gsl_vector* val;
  gsl_matrix* vec;
  gsl_eigen_symmv_workspace* w;
  int i, status;

  val = gsl_vector_alloc(3);
  vec = gsl_matrix_alloc(3, 3);
  w = gsl_eigen_symmv_alloc(3);

  i = 0;
  B[i++] = A[XX];
  B[i++] = A[XY];
  B[i++] = A[XZ];
  B[i++] = A[YX];
  B[i++] = A[YY];
  B[i++] = A[YZ];
  B[i++] = A[ZX];
  B[i++] = A[ZY];
  B[i++] = A[ZZ];

  m = gsl_matrix_view_array(B, 3, 3);
  status = gsl_eigen_symmv(&m.matrix, val, vec, w);
  if (status != GSL_SUCCESS) ERR(("gsl_eigen_symmv failed"));
  gsl_eigen_symmv_sort(val, vec, GSL_EIGEN_SORT_ABS_ASC);

  i = 0;
  VAL[i] = gsl_vector_get(val, i);
  i++;
  VAL[i] = gsl_vector_get(val, i);
  i++;
  VAL[i] = gsl_vector_get(val, i);

  gsl_vector_free(val);
  gsl_matrix_free(vec);
  gsl_eigen_symmv_free(w);
  return 0;
}

int chmath_eig_vectors(
    const double A[6], /**/ double a[3], double b[3], double c[3]) {
  enum { XX, XY, XZ, YY, YZ, ZZ };
  enum { YX = XY, ZX = XZ, ZY = YZ };
  enum { X, Y, Z };

  double B[3 * 3];
  gsl_matrix_view m;
  gsl_vector* val;
  gsl_matrix* vec;
  gsl_eigen_symmv_workspace* w;
  int i, status;

  val = gsl_vector_alloc(3);
  vec = gsl_matrix_alloc(3, 3);
  w = gsl_eigen_symmv_alloc(3);

  i = 0;
  B[i++] = A[XX];
  B[i++] = A[XY];
  B[i++] = A[XZ];
  B[i++] = A[YX];
  B[i++] = A[YY];
  B[i++] = A[YZ];
  B[i++] = A[ZX];
  B[i++] = A[ZY];
  B[i++] = A[ZZ];

  m = gsl_matrix_view_array(B, 3, 3);
  status = gsl_eigen_symmv(&m.matrix, val, vec, w);
  if (status != GSL_SUCCESS) ERR(("gsl_eigen_symmv failed"));
  gsl_eigen_symmv_sort(val, vec, GSL_EIGEN_SORT_ABS_ASC);

  i = 0;
  a[X] = get(vec, X, X);
  a[Y] = get(vec, Y, X);
  a[Z] = get(vec, Z, X);
  b[X] = get(vec, X, Y);
  b[Y] = get(vec, Y, Y);
  b[Z] = get(vec, Z, Y);
  c[X] = get(vec, X, Z);
  c[Y] = get(vec, Y, Z);
  c[Z] = get(vec, Z, Z);
  gsl_vector_free(val);
  gsl_matrix_free(vec);
  gsl_eigen_symmv_free(w);
  return 0;
}
