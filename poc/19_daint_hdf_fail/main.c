#include <mpi.h>
#include <hdf5.h>
#include <stdio.h>

#define WARN(x) do {					\
	fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
	dprint x;					\
	fputs("\n", stderr);				\
  } while (0)
static int
dprint(const char *fmt, ...)
{
  int r;
  va_list ap;

  va_start(ap, fmt);
  r = vfprintf(stderr, fmt, ap);
  va_end(ap);
  return r;
}


int main() {
  unsigned a, b, c;
  H5get_libversion(&a, &b, &c);

  hid_t plist, file;
  char full[100];

  plist = H5Pcreate(H5P_FILE_ACCESS);
  if (plist < 0) {
    WARN(("H5Pcreate failed"));
    return plist;
  }
  H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
}
