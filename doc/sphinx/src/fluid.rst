.. _s:fluid:

Navier-Stokes equations
=======================


Abstract class ``FluidSolver`` describes the interface
of a solver for the Navier-Stokes equations

.. math::
  \nabla \cdot \mathbf{u} &= S_v
  \\
  \rho \Big(
  \frac{\partial \mathbf{u}}{\partial t}
  + (\mathbf{v}\cdot\nabla) \mathbf{u}
  \Big)
  &=  -\nabla p + \nabla \cdot \mu (\nabla \mathbf{u} + \nabla \mathbf{u}^T)
  + \mathbf{f}

where the continuity equation is discretized as

.. math::
  \sum_{f\in c} v_f S^c_{\!f} = S_v

and the momentum equation as :eq:`e:momentum`,
where :math:`S_v` is the volume source term.

The constructor takes pointers to fields
of density :math:`\rho`,
viscosity :math:`\mu`,
force :math:`\mathbf{f}`
and volume source :math:`S_v`:

.. includecode:: src/solver/fluid.h
  :func: FluidSolver
  :comment:
  :dedent: 2

After an implementation is constructed, the solution
is advanced by calling ``MakeIteration()``
and the current solution is provided by

.. includecode:: src/solver/fluid.h
  :func: GetVelocity
  :dedent: 2

.. includecode:: src/solver/fluid.h
  :func: GetPressure
  :dedent: 2

.. includecode:: src/solver/fluid.h
  :func: GetVolumeFlux
  :dedent: 2

The volume flux satisfies the continuity equation
and the cell-based velocity satisfies the momentum equation.

SIMPLE
------

Class ``Simple`` implements the solver using the SIMPLE method.

In addition to the requirements of the base class,
the constructor takes mappings describing the boundary
conditions and the initial fields.

.. includecode:: src/solver/simple.h
  :func: Simple
  :dedent: 2

with the initial velocity ``fcw``,
face conditions ``mfc``,
cell conditions ``mcc``,
density ``fcr``,
viscosity ``fcd``,
force ``fcf``,
projections of well-balanced force ``fcbp``,
volume source ``fcsv`` and mass source ``fcsm``.
Parameters of the solver are provided by

.. includecode:: src/solver/simple.h
  :struct: SimplePar
  :dedent: 0

Projection
----------

Class ``Proj`` implements the solver using Chorin's projection.

In addition to the requirements of the base class,
the constructor takes mappings describing the boundary
conditions and the initial fields.

.. includecode:: src/solver/proj.h
  :func: Proj
  :dedent: 2

with the initial velocity ``fcw``,
face conditions ``mfc``,
cell conditions ``mcc``,
density ``fcr``,
viscosity ``fcd``,
force ``fcf``,
projections of well-balanced force ``fcbp``,
volume source ``fcsv`` and mass source ``fcsm``.
Parameters of the solver are provided by

.. includecode:: src/solver/proj.h
  :struct: ProjPar
  :dedent: 0

Boundary conditions
-------------------

The boundary conditions are specified by a map
from ``IdxFace`` to ``CondFaceFluid``,
instances of which can be generated with function ``solver::Parse()``
from a string.
Value of ``id`` is written to field ``fluidcond`` in ``bc.vtk``
(see :ref:`s:output`).

.. table:: Fluid boundary conditions.
   :name: t:fluid_boundary

   +---------------------+--------------------------+-----------------------------------+----+
   | class               | Parse() format           | description                       | id |
   +=====================+==========================+===================================+====+
   | ``NoSlipWallFixed`` | ``wall <x y z>``         | no-slip wall with fixed velocity  |  1 |
   +---------------------+--------------------------+-----------------------------------+----+
   | ``InletFixed``      | ``inlet <x y z>``        | inlet with given velocity         |  3 |
   +---------------------+--------------------------+-----------------------------------+----+
   | ``InletFlux``       | ``inletflux <x y z id>`` | inlet with given total flux       |  3 |
   +---------------------+--------------------------+-----------------------------------+----+
   | ``OutletAuto``      | ``outlet``               | outlet                            |  4 |
   +---------------------+--------------------------+-----------------------------------+----+
   | ``SlipWall``        | ``slipwall``             | free-slip wall                    |  2 |
   +---------------------+--------------------------+-----------------------------------+----+
