#include <stdio.h>
#include <hdf5.h>

#ifndef H5_HAVE_PARALLEL
#error needs parallel HDF5
#endif

int
main(void)
{
  	unsigned int maj, min, rel;
	H5get_libversion(&maj, &min, &rel);
        fprintf(stderr, "hdf5: %d.%d.%d\n", maj, min, rel);
}
