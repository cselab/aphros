.. marching cubes

Marching cubes
==============

SYNOPSIS
--------

Marching cubes.

::

   #include <marh.h>

   static double MARCH_O[][3] = {
     {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
     {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
   };
   enum { MARCH_NTRI = 48 };
   int march_cube(double u[8], int* ntri, double* tri);
   int march_cube_location(double u[8], int* ntri, double* tri, int* first, int* second, double* offset);

DESCRIPTION
-----------

Generates triangles of the isosurface for three-dimensional discrete
field.

march_cube generates triangles for one cell, u[8] are values of the
field in the nodes. Returns ntri -- the number of triangles, tri --
vertices of the triangles in the form [x0 y0 z0 ... xn yn zn].

::

   double cube[8] = {-1, 0, 0, 0, 0, 0, 0, 1};
   double tri[3 * 3 * MARCH_NTRI];
   int n;

   march_cube(cube, &n, tri);

MARCH_NTRI is a maximum number of triangles.
MARCH_O[8][3] are coordinates of the nodes.

march_cube_location in addition to the output of the march_cube
returns location of the vertex i relative the nodes.
first[i], second[i] are indexes of the nodes. offset[i] is a
relative distance from the node first[i].

::

   enum {X, Y, Z};
   static double av(double a, double b, double o) {
     return a + (b - a) * o;
   }

   double cube[8] = {-1, 0, 0, 0, 0, 0, 0, 1};
   double tri[3 * 3 * MARCH_NTRI];
   int first[3 * MARCH_NTRI], second[3 * MARCH_NTRI];
   double offset[3 * MARCH_NTRI];
   double x;
   int u;

   u = 0;
   march_cube_location(cube, &n, tri, first, second, offset);
   x = av(MARCH_O[first[u]][X], MARCH_O[second[u]][X], offset[u]);
   assert(fabs(pos - tri[3 * u + X] < 1e-12));

SOURCE
------

.. linkpath:: deploy/lib/march

EXAMPLES
--------

.. linkpath:: deploy/lib/march/example/sphere.c
.. linkpath:: deploy/lib/march/example/main.c
