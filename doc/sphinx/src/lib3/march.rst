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

Use marching cube algorithm to generate triangles of the isosurface
for three-dimensional discrete field.

march_cube generate triangles for one cubic cell, u[8] are values of
the field in the corners of the cell. Returns ntri -- the number of
triangles, tri -- vertices of the triangles in the form [x0 y0 z0
... xn yn zn].

::

   double cube[8] = {-1, 0, 0, 0, 0, 0, 0, 1};
   double tri[3 * 3 * MARCH_NTRI];
   int n;

   march_cube(cube, &n, tri);

MARCH_NTRI is a maximum number of triangles.
MARCH_O[8][3] are coordinates of the corners.

march_cube_location in addition to the output of the march_cube
returns location of the vertices i relative the corners of the cell.
first[i], second[i] are indexes of the cell's corners. offset[i] is a
relative distance from the corner first[i].

::

   enum {X, Y, Z};
   static double av(double a, double b, double o) {
     return a + (b - a) * o;
   }

   double cube[8] = {-1, 0, 0, 0, 0, 0, 0, 1};
   double tri[3 * 3 * MARCH_NTRI];
   int p[3 * MARCH_NTRI], q[3 * MARCH_NTRI];
   double offset[3 * MARCH_NTRI];
   double pos;
   int u;

   u = 0;
   march_cube_location(cube, &n, tri, p, q, offset);
   pos = av(MARCH_O[p[u]][X], MARCH_O[q[u]][X], offset[u]);
   assert(fabs(pos - tri[3 * u + X] < 1e-12));
