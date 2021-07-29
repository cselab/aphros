.. _s:output:

Output
======

Simulations produce output of various types:

* statistics (stdout)
* fields (HDF5, plain or XML VTK)
* PLIC polygons (polydata legacy VTK)
* marching cubes triangles (polydata legacy VTK)
* boundary conditions (polydata legacy VTK)

Boundary conditions are written to ``bc.vtk`` if enabled by ``set int
dumpbc 1``.  It contains cell fields ``block`` (block id), ``cond``
(conditions for advection) and ``condfluid`` (conditions for fluid,
see :ref:`t:fluid_boundary`).

The output is implemented in function

.. includecode:: src/util/hydro.ipp
  :func: DumpBcPoly

