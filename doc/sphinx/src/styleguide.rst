Style guide
===========

Use `Google C++ Style Guide <https://google.github.io/styleguide/cppguide.html>`_.

Naming conventions
------------------

-  Variable:

   -  lower case with underscores
   -  single letter when possible
   -  if first letter, put that word in comment: n // name
   -  if another letter, put that word with letter in […]: a // n[a]me
   -  up to 4 letters if possible
   -  private variables end with underscore: ``a_``; exceptions: m
      (Mesh), var (Vars)

-  Template argument:

   -  CamelCase for types and lower case for values
   -  if alias needed, ends with underscore: ``Mesh_``
   -  else one letter

-  Filename:

   -  lower case, extensions h and cpp

-  Class

   -  …

-  No names starting with underscore, reserved
-  For indices: ``i`` - generic ``c`` - IdxCell ``f`` - IdxFace ``w`` -
   MIdx ``d`` - direction ``d`` - index 0..dim
-  Other: ``u`` - generic field ``q`` - neighbour cell id ``o`` - other
-  Treat pointers same as values


Toolbox classes
---------------

-  combines multiple templated functions to avoid repeating template
   parameters
-  acts as namespace
-  contains only static functions

Comments
--------

-  capitalized descriptive for a class or function (e.g. Creates
   instance)
-  capitalized imperative for expressions in implementation (e.g. Create
   instance)
-  non-capitalized for declarations (e.g. buffer index)

Units
-----

Example: ``Scal s; // source [density/time]``

-  length
-  time
-  mass
-  velocity (length/time)
-  area (length^2)
-  volume (length^3)
-  density (mass/volume)
-  value (user defined, e.g. 1 for concentration)

Ordering
~~~~~~~~

-  parent ``P``
-  kernel ``K``
-  mesh ``M``
-  constructor parameters ``Par``
-  dimension ``dim``
