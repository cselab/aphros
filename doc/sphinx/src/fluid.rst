.. _s:fluid:

Navier-Stokes equations
=======================


Abstract class ``FluidSolver`` describes the interface
of a solver for the Navier-Stokes equations

.. math::
  \nabla \cdot \mathbf{u} = S_v
  \\
  \rho \Big(
  \frac{\partial \mathbf{u}}{\partial t}
  + (\mathbf{v}\cdot\nabla) \mathbf{u}
  \Big)
  + \nabla \cdot (\mu \nabla \mathbf{u})
  + \mathbf{f}
  = 0

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
  :struct: Par
