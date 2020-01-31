.. _s:primlist:

Geometric primitives
====================

Header ``func/primlist.h`` provides routines to define level-set functions
for a list of geometric primitives such as spheres and boxes.
A primitive is described by

.. includecode:: src/func/primlist.h
  :struct: Primitive

Function ``Parse()`` parses a stream with a list of primitives.

.. includecode:: src/func/primlist.h
  :func: Parse
  :comment:

Available primitives and their parameters:

``sphere``
  ``cx cy cz`` (center), ``rx ry rz`` (half-size);
``box``
  ``cx cy cz`` (center), ``rx ry rz`` (half-size);
``ring``
  ``cx cy cz`` (center), ``nx ny nz`` (normal), ``r`` (radius), ``th`` (thickness);
``smooth_step``
  ``cx cy cz`` (center), ``nx ny nz`` (normal), ``tx ty tz`` (tangent),
  ``ln`` (size along normal), ``lt`` (size along tangent);
``cylinder``
  ``cx cy cz`` (center), ``tx ty tz`` (axis), ``r`` (radius),
  ``t0 t1`` (coordinates relative to center).

Each primitive defines a level-set function which is positive inside the body.
By default, the resulting level-set function is composed from the list of
primitives using the union operation (taking the maximum value).
To change the default operation, modifiers can be added
before the name of the primitive:

* ``-``: minus, multiply level-set by -1;
* ``&``: intersection, take minimum with level-set from all entries above.

Example of a list of primitives

.. includecode:: examples/200_primlist/b.dat
  :language: none

(see setup in :linkpath:`examples/200_primlist`)

.. image:: ../../../examples/200_primlist/a.jpg
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
