.. marching cubes

Marching cubes
==============

Synopsis
--------

Marching cubes.

.. code-block:: c++

   #include <march.h>

   static double MARCH_O[][3] = {
     {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
     {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
   };
   enum { MARCH_NTRI = 48 };
   int march_cube(double u[8], int* ntri, double* tri);
   int march_cube_location(double u[8], int* ntri, double* tri, int* first, int* second, double* offset);

Description
-----------

Generates triangles of the isosurface for three-dimensional discrete
field.

march_cube generates triangles for one cell, u[8] are values of the
field in the nodes. Returns ntri -- the number of triangles, tri --
vertices of the triangles in the form [x0 y0 z0 ... xn yn zn].

.. code-block:: c++

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

.. includecode:: examples/000_march/main.c

Source
------

:linkpath:`deploy/lib/march`

Examples
--------

Generate triangles for cube configurations given as arguments:

| :linkpath:`examples/001_march_cell/main.c`

Generate a triangulated sphere:

| :linkpath:`examples/002_march_sphere/main.c`
