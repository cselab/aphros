/*! gcc -Wall -g -o test test.c libkdtree.a */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include "kdtree.h"

unsigned int get_msec(void)
{
	static struct timeval timeval, first_timeval;

	gettimeofday(&timeval, 0);

	if(first_timeval.tv_sec == 0) {
		first_timeval = timeval;
		return 0;
	}
	return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}


int main(int argc, char **argv)
{
	int i, vcount = 10;
	void *kd, *set;
	unsigned int msec, start;

	if(argc > 1 && isdigit(argv[1][0])) {
		vcount = atoi(argv[1]);
	}
	printf("inserting %d random vectors... ", vcount);
	fflush(stdout);

	kd = kd_create(3);

	start = get_msec();
	for(i=0; i<vcount; i++) {
		float x, y, z;
		x = ((float)rand() / RAND_MAX) * 200.0 - 100.0;
		y = ((float)rand() / RAND_MAX) * 200.0 - 100.0;
		z = ((float)rand() / RAND_MAX) * 200.0 - 100.0;

		assert(kd_insert3(kd, x, y, z, 0) == 0);
	}
	msec = get_msec() - start;
	printf("%.3f sec\n", (float)msec / 1000.0);

	start = get_msec();
	set = kd_nearest_range3(kd, 0, 0, 0, 40);
	msec = get_msec() - start;
	printf("range query returned %d items in %.5f sec\n", kd_res_size(set), (float)msec / 1000.0);
	kd_res_free(set);

	kd_free(kd);
	return 0;
}
