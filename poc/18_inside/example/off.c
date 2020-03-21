#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <inside.h>

const char* me = "off";

static void usg(void) {
  fprintf(stderr, "%s < off > xyz\n", me);
  exit(1);
}

int
main(int argc, const char** argv) {
  int nt;
  int nv;
  int *tri;
  double *x;
  double *y;
  double *z;
  int status;

  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    default:
      fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
      exit(2);
    }
  off_read(stdin, &status, &nt, &tri, &nv, &x, &y, &z);
  off_write(nt, tri, nv, x, y, z, stdout);
  off_fin(tri, x, y, z);
}

enum {SIZE = 999};
int
off_read(FILE * f, int *status, int * pnt, int ** ptri, int * pnv, double ** px, double ** py, double ** pz)
{
  enum {Off, Numbers, Ver, Tri, End};
  int state;
  char line[SIZE];
  char *s;
  int i;
  int j;
  int v;
  int t;
  int nv;
  int nt;
  int *tri;
  int npt;
  int t0;
  int t1;
  int t2;
  int cnt;
  double *x;
  double *y;
  double *z;

  state = Off;
  v = t = 0;
  for (;;) {
    if ((s = fgets(line, SIZE, f)) == NULL)
      break;
    while (isspace(*s)) s++; /* leading spaces */
    for (i = 0; s[i] != '\0'; i++)
      if (s[i] == '\n' || s[i] == '#') {
	s[i] = '\0';
	break;
      }
    for (i = j = 0; s[i] != '\0'; i++) /* trailing spaces */
      if (!isspace(s[i]))
	j = i;
    s[j + 1] = '\0';
    if (s[0] == '\0') /* empty line */
      continue;
    switch (state) {
    case Off:
      if (strncmp(s, "OFF", SIZE))
	goto not_off;
      state = Numbers;
      break;
    case Numbers:
      if (sscanf(s, "%d %d %*d", &nv, &nt) != 2) {
	fprintf(stderr, "%s:%d: expecting numbers, got '%s'\n", s, __FILE__, __LINE__);
	goto err;
      }
      x = malloc(nv*sizeof(*x));
      y = malloc(nv*sizeof(*y));
      z = malloc(nv*sizeof(*z));
      tri = malloc(3*nt*sizeof(*tri));
      if (tri == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	goto err;
      }
      state = Ver;
      break;
    case Ver:
      if (sscanf(s, "%f %f %f", &x[v], &y[v], &z[v]) != 3) {
	fprintf(stderr, "%s:%d: expcting vertices got '%s'\n", __FILE__, __LINE__, s);
	goto err;
      }
        fprintf(stderr, "%g %g %g\n", x[v], y[v], z[v]);
      exit(2);
      v++;
      if (v == nv)
	state = Tri;
      break;
    case Tri:
      cnt = sscanf(s, "%d %d %d %d", &npt, &t0, &t1, &t2);
      if (cnt != 4 || npt != 3) {
	fprintf(stderr, "%s:%d: expcting triangle got '%s'\n", __FILE__, __LINE__, s);
	goto err;
      }
      tri[t++] = t0;
      tri[t++] = t1;
      tri[t++] = t2;
      if (t == 3*nt)
	state = End;
      break;
    case End:
      fprintf(stderr, "%s:%d: extra line '%s'\n", __FILE__, __LINE__, s);
      goto err;
      break;
    }
  }
  if (state != End) {
      fprintf(stderr, "%s:%d: off file is not complite\n", __FILE__, __LINE__);
      goto err;
  }
  *pnv = nv;
  *pnt = nt;
  *ptri = tri;
  *px = x;
  *py = y;
  *pz = z;
  *status = 0;
  return 0;
 err:
  return 1;
 not_off:
  *status = 0;
  return 0;
}

int
off_fin(int *tri, double *x, double *y, double *z)
{
  free(tri);
  free(x);
  free(y);
  free(z);
}

int
off_write(int nt, int * tri, int nv, double * x, double * y, double * z, FILE *f)
{
  int i;
  if (fprintf(f, "OFF\n") < 0) {
      fprintf(stderr, "%s:%d: fail to write\n", __FILE__, __LINE__);
      goto err;
  }
  fprintf(f, "%d %d 0\n", nv, nt);
  for (i = 0; i < nv; i++)
    fprintf(f, "%.16g %.16g %.16g\n", x[i], y[i], z[i]);
  for (i = 0; i < nt; i++)
    fprintf(f, "3 %d %d %d\n", tri[3*i], tri[3*i+1], tri[3*i+2]);
 err:
  return 1;
}
