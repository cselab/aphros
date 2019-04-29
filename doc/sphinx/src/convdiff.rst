Convection-diffusion equation
=============================

Abstract class ``ConvDiffVect`` describes the interface
of a solver for the vector convection-diffusion equation

.. math::
  \rho \Big(
  \frac{\partial \mathbf{u}}{\partial t} 
  + (\mathbf{v}\cdot\nabla) \mathbf{u}
  \Big)
  + \nabla \cdot (\mu \nabla \mathbf{u})
  + \mathbf{f}
  = 0

in the discrete form

.. math::
  \rho_c \Big(
  \frac{\delta \mathbf{u}_c}{\delta t}
  + \frac{1}{V_c}\sum_{f\in c} v_f \mathbf{u}_f S^c_{\!f}
  \Big)
  + \frac{1}{V_c}\sum_{f\in c} \mu_f
    \frac{\delta \mathbf{u}_f}{\delta n} S^c_{\!f}
  + \mathbf{f}_c = 0,

where :math:`c` is a cell index,
:math:`\phi_c` is a cell-average,
:math:`f\in c` are neighbour faces of cell :math:`c`,
:math:`\phi_f` is a face-average,
:math:`S^c_{\!f}` is the signed face area
(positive if :math:`\mathbf{n}_f` is an outer normal to :math:`c`)
and :math:`V_c` is the cell volume.

.. .. literalinclude:: src/solver/convdiffi.ipp
  :language: cpp
  :lines: 61-136

The constructor takes pointers to fields
of density :math:`\rho`,
viscosity :math:`\mu`,
force :math:`\mathbf{f}`
and volume flux :math:`\mathbf{v} \cdot \mathbf{S}_f`.

.. includecode:: src/solver/convdiffv.h
  :func: ConvDiffVect
  :comment:
  :dedent: 2

After an instance is constructed, the solution 
is advanced by calling the inherited ``MakeIteration()``
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
  :func: GetVelocityEquations 
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
  + \frac{1}{V_c}\sum_{f\in c} \mu_f^s
    \frac{\delta \mathbf{u}_f^{s+1}}{\delta n} S^c_{\!f}
  + \mathbf{f}_c^s = 0.

which requires solving a linear system at every iteration.
The system is solved by 

.. includecode:: src/solver/convdiffi.ipp
  :func: Solve
  :comment:
  :dedent: 2

creating a request for 

.. includecode:: src/geom/mesh.h
  :func: Solve
  :dedent: 2

which leads to calling the Hypre library

.. includecode:: src/linear/hypre.h
  :func: Hypre
  :dedent: 2


Explicit solver
---------------

Class ``ConvDiffVectExp`` implements an explicit solver
corresponding to the discrete equation

.. math::
  \rho_c^s \Big(
  \frac{\delta \mathbf{u}_c^{s+1}}{\delta t}
  + \frac{1}{V_c}\sum_{f\in c} v_f^s \mathbf{u}_f^{s} S^c_{\!f}
  \Big)
  + \frac{1}{V_c}\sum_{f\in c} \mu_f^s
    \frac{\delta \mathbf{u}_f^{s}}{\delta n} S^c_{\!f}
  + \mathbf{f}_c^s = 0.

Here the solution is advanced by explicit formulas
and the linear system is constructed
only to implement ``Assemble()`` and 
``GetVelocityEquations()`` of the interface.
