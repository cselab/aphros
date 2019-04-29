Convection-diffusion solver
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

The following function returns component `d` of the vector equations
in terms of velocity correction :math:`u'=u^{s+1}-u^s`:

.. includecode:: src/solver/convdiffv.h
  :func: GetVelocityEquations 
  :dedent: 2

