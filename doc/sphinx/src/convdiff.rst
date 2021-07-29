.. _s:convdiff:

Convection-diffusion equation
=============================


Abstract class ``ConvDiffVect`` describes the interface
of a solver for the vector convection-diffusion equation

.. math::
  \rho \Big(
  \frac{\partial \mathbf{u}}{\partial t}
  + (\mathbf{v}\cdot\nabla) \mathbf{u}
  \Big)
  = \nabla \cdot (\mu \nabla \mathbf{u})
  + \mathbf{f}

in the discrete form

.. math::
  :label: e:momentum

  \rho_c \Big(
  \frac{\delta \mathbf{u}_c}{\delta t}
  + \frac{1}{V_c}\sum_{f\in c} v_f \mathbf{u}_f S^c_{\!f}
  \Big)
  = \frac{1}{V_c}\sum_{f\in c} \mu_f
    \frac{\delta \mathbf{u}_f}{\delta n} S^c_{\!f}
  + \mathbf{f}_c

where :math:`c` is a cell index,
:math:`\phi_c` is a cell-average,
:math:`f\in c` are neighbour faces of cell :math:`c`,
:math:`\phi_f` is a face-average,
:math:`S^c_{\!f}` is the signed face area
(positive if :math:`\mathbf{n}_f` is an outer normal to :math:`c`)
and :math:`V_c` is the cell volume.

The constructor takes pointers to fields
of density :math:`\rho`,
viscosity :math:`\mu`,
force :math:`\mathbf{f}`
and volume flux :math:`\mathbf{v} \cdot \mathbf{S}_f`:

.. includecode:: src/solver/convdiffv.h
  :func: ConvDiffVect
  :comment:
  :dedent: 2

After an implementation is constructed, the solution
is advanced by calling ``MakeIteration()``
and the current velocity field is returned by

.. includecode:: src/solver/convdiffv.h
  :func: GetVelocity
  :comment:
  :dedent: 2

The interface also exposes functions

.. includecode:: src/solver/convdiffv.h
  :func: Assemble
  :dedent: 2

.. includecode:: src/solver/convdiffv.h
  :func: GetDiag
  :dedent: 2

.. includecode:: src/solver/convdiffv.h
  :func: GetConst
  :dedent: 2

to assemble the linear system for a given velocity field
``fcw`` from the previous iteration and the volume flux ``ffv``
and access its scalar component ``d``.
Furthermore, the current solution can be explicitly corrected
with

.. includecode:: src/solver/convdiffv.h
  :func: CorrectVelocity
  :dedent: 2

Both these features are required for the pressure correction equation.


Implicit solver
---------------

Class ``ConvDiffVectImp`` implements an implicit solver
corresponding to the discrete equation

.. math::
  \rho_c^s \Big(
  \frac{\delta \mathbf{u}_c^{s+1}}{\delta t}
  + \frac{1}{V_c}\sum_{f\in c} v_f^s \mathbf{u}_f^{s+1} S^c_{\!f}
  \Big)
  = \frac{1}{V_c}\sum_{f\in c} \mu_f^s
    \frac{\delta \mathbf{u}_f^{s+1}}{\delta n} S^c_{\!f}
  + \mathbf{f}_c^s

which requires solving a linear system at every iteration.

In addition to the requirements of the base class,
the constructor takes mappings describing the boundary
conditions and the initial fields:

.. includecode:: src/solver/convdiffv.h
  :func: ConvDiffVect
  :dedent: 2

with the initial velocity ``fcvel``,
face conditions ``mfc``,
cell conditions ``mcc``,
density ``fcr``,
viscosity ``fcd``,
force ``fcs``,
and volume flux ``ffv``.
Parameters of the solver are provided by

.. includecode:: src/solver/convdiff.h
  :struct: ConvDiffPar
  :dedent: 0

where ``sc`` defines the interpolation scheme

.. table:: Interpolation schemes.
  :name: t:interp_schemes

   +--------------------+---------------------------------+
   | ``ConvSc::fou``    | First Order Upwind              |
   +--------------------+---------------------------------+
   | ``ConvSc::cd``     | Central Differences (midpoint)  |
   +--------------------+---------------------------------+
   | ``ConvSc::sou``    | Second Order Upwind             |
   +--------------------+---------------------------------+
   | ``ConvSc::quick``  | QUICK                           |
   +--------------------+---------------------------------+




Explicit solver
---------------

Class ``ConvDiffVectExp`` implements an explicit solver
corresponding to the discrete equation

.. math::
  \rho_c^s \Big(
  \frac{\delta \mathbf{u}_c^{s+1}}{\delta t}
  + \frac{1}{V_c}\sum_{f\in c} v_f^s \mathbf{u}_f^{s} S^c_{\!f}
  \Big)
  = \frac{1}{V_c}\sum_{f\in c} \mu_f^s
    \frac{\delta \mathbf{u}_f^{s}}{\delta n} S^c_{\!f}
  + \mathbf{f}_c^s.

Here the solution is advanced by explicit formulas
and the linear system is constructed
only to implement ``Assemble()`` and
``GetVelocityEquations()`` of the base class.
