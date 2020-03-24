/*! gcc -std=c89 -pedantic -Wall -g -o test2 test2.c libkdtree.a -lm */
/* Extended test program, contributed by David Underhill */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include "kdtree.h"

#define DEF_NUM_PTS 10

/* returns the distance squared between two dims-dimensional double arrays */
static double dist_sq( double *a1, double *a2, int dims );

/* get a random double between -10 and 10 */
static double rd( void );

int main(int argc, char **argv) {
  int i, num_pts = DEF_NUM_PTS;
  void *ptree;
  int *data, pch;
  struct kdres *presults;
  double pos[3], dist;
  double pt[3] = { 0, 0, 1 };
  double radius = 10;

  if(argc > 1 && isdigit(argv[1][0])) {
    num_pts = atoi(argv[1]);
  }

  if(!(data = malloc(num_pts*sizeof(*data)))) {
    perror("malloc failed");
    return 1;
  }

  srand( time(0) );
  ptree = kd_create();

  /* add some random nodes to the tree (assert nodes are successfully inserted) */
  double r[3];
  for( i=0; i<num_pts; i++ ) {
    data[i] = i;
    r[0] = rd();
    r[1] = rd();
    r[2] = rd();
    assert( 0 == kd_insert( ptree, r, data[i] ) );
  }

  /* find points closest to the origin and within distance radius */
  presults = kd_nearest_range( ptree, pt, radius );

  /* print out all the points found in results */
  printf( "found %d results:\n", kd_res_size(presults) );

  while( !kd_res_end( presults ) ) {
    /* get the data and position of the current result item */
    pch = kd_res_item( presults, pos );

    /* compute the distance of the current result from the pt */
    dist = sqrt( dist_sq( pt, pos, 3 ) );

    /* print out the retrieved data */
    printf( "node at (%.3f, %.3f, %.3f) is %.3f away and has data=%d\n",
	    pos[0], pos[1], pos[2], dist, pch );

    /* go to the next entry */
    kd_res_next( presults );
  }

  /* free our tree, results set, and other allocated memory */
  free( data );
  kd_res_free( presults );
  kd_free( ptree );

  return 0;
}

static double dist_sq( double *a1, double *a2, int dims ) {
  double dist_sq = 0, diff;
  while( --dims >= 0 ) {
    diff = (a1[dims] - a2[dims]);
    dist_sq += diff*diff;
  }
  return dist_sq;
}

static double rd( void ) {
  return (double)rand()/RAND_MAX * 20.0 - 10.0;
}