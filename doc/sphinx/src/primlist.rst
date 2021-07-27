.. _s:primlist:

Geometric primitives
====================

Header ``func/primlist.h`` provides routines to define level-set functions
for a list of geometric primitives such as spheres and boxes.
A primitive is described by

.. includecode:: src/func/primlist.h
  :struct: Primitive

Function ``GetPrimitives()`` parses a stream with a list of primitives.

.. includecode:: src/func/primlist.h
  :func: GetPrimitives
  :comment:

Available primitives and their parameters:

``sphere``
  Ellipsoid with principal axes aligned with the coordinate axes.
  Parameters:
  ``cx cy cz`` (center), ``rx ry rz`` (half-size);
``box``
  Rectangular box with sides aligned with the coordinate planes.
  Parameters:
  ``cx cy cz`` (center), ``rx ry rz`` (half-size);
``ring``
  Torus.
  Parameters:
  ``cx cy cz`` (center), ``nx ny nz`` (normal), ``r`` (radius), ``th`` (thickness);
``smooth_step``
  Smooth step :cite:t:`almgren1997cartesian`.
  Parameters:
  ``cx cy cz`` (center), ``nx ny nz`` (normal), ``tx ty tz`` (tangent),
  ``ln`` (size along normal), ``lt`` (size along tangent);
``cylinder``
  Right circular cylinder.
  Parameters:
  ``cx cy cz`` (center), ``nx ny nz`` (normal), ``r`` (radius),
  ``n0 n1`` (range along normal relative to center);
``polygon``
  Cylinder bounded by parallel planes with the plane section specified as a
  sequence of non-intersecting polygons.
  Parameters:
  ``ox oy oz`` (origin), ``nx ny nz`` (normal), ``ux uy uz``
  (direction of 2D x-axis), ``n0 n1`` (range along normal relative to origin),
  ``scale`` (factor applied to 2D vertices), ``x y ...``
  (2D vertices of all polygons, first and last vertices of each polygon must
  coincide);

``ruled``
  Ruled surface bounded by parallel planes with two plane sections
  on the opposite sides specified as two sequences of non-intersecting polygons.
  Parameters:
  ``ox oy oz`` (origin), ``nx ny nz`` (normal), ``ux uy uz`` (direction of 2D
  x-axis), ``n0 n1`` (range along normal relative to origin, sides 0 and 1),
  ``scale0 scale1``
  (factors applied to 2D vertices on sides 0 and 1), ``x y ...``
  (2D vertices of all polygons on side 0, first and last vertices of
  each polygon must coincide), ``x y ...`` (same but on side 1).

Each primitive defines a level-set function which is positive inside the body.
By default, the resulting level-set function is composed from the list of
primitives using the union operation (taking the maximum value).
To change the default operation, modifiers can be added
before the name of the primitive:

* ``-``: minus, multiply level-set by -1;
* ``&``: intersection, take the minimum with the current level-set.

Example of a list of primitives

.. includecode:: examples/200_primlist/b.dat
  :language: none

(see setup in :linkpath:`examples/200_primlist`)

.. image:: ../../../examples/200_primlist/a_orig.jpg
  :width: 400
  :align: center

Parameter ``list_ls`` in the configuration defines the
method of computing the volume fraction field from the level-set functions

``set int list_ls 0``
  step-wise approximation (1 if level-set is positive, 0 otherwise);
``set int list_ls 1``
  linear approximation with normal and plane constant
  computed from the level-set at the cell center,
  does not support modifiers;
``set int list_ls 2``
  using the ``overlap`` library to compute the exact
  volume fraction cut by an ellipsoid (only valid for primitive ``sphere``);
``set int list_ls 3``
  linear approximation with normal and plane constant
  computed from the level-set on mesh nodes, supports modifiers.
