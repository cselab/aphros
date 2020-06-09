.. _s:vofm:

Multilayer VOF
==============

Multilayer VOF is a method for solving the advection equation
that can represent multiple interfaces in the same cell.
This enables simulations of foams, suspensions,
and other multiphase systems where bubbles and drops that do not coalesce.

The split conservative advection scheme [weymouth2010]_ operates on a volume
fraction field and uses the PLIC reconstruction to compute the fluxes.
In three dimensions, one advection step is split into three sub-steps, or
sweeps. To ensure conservation, a divergence term is added explicitly.
The normals are estimated using the mixed Youngs-centered scheme.

The multilayer VOF implements the same advection scheme
but operates on *multilayer fields*.
A multilayer field *(u,q)* stores mappings between *colors* *q* and *values* *u*.
The colors are real numbers, and *none=-1* denotes an empty color.
Assuming that each cell can contain at most *L=4* values,
the mappings are represented as *L*-tuples up to reordering.
For example, the following tuples encode a field *(u,q)*
with three values in cell *c1* and two values in cell *c2*

::

  u[c1]: {0.2, 0.4,  0.1,    0}
  q[c1]: {  0,   1,    2, none}
  u[c2]: {0.3, 0.7,    0,    0}
  q[c2]: {  1,   2, none, none}

The ordering of tuples is insignificant, so simultaneous permutation
of values and colors gives equivalent representations

::

  u[c1]: {0.4, 0.1, 0.2,    0}
  q[c1]: {  1,   2,   0, none}
  u[c2]: {0.7,    0, 0.3,    0}
  q[c2]: {  2, none,   1, none}

A conventional field can be constructed from a multilayer field
by selecting values with a given color.
For example, selecting the values with color *qs=2* gives a conventional field
*us*

::

  us[c1]: 0
  us[c2]: 0.7

The following pseudocode implements the multilayer VOF method.
One advection step is split in the spatial directions
and in three dimensions consists of three substeps, or sweeps.

.. code-block:: c++

    void Step() {
      for (d : {0, 1, 2}) {
        {u, q} = Sweep(d, u, v, q, n, a);
        {n, a} = Plic(u);
      }
      q = Recolor(u, q);
    }

Given the PLIC reconstruction of the interface,
one sweep computes fluxes, updates the volume fractions
and propagates the colors to downwind cells

.. code-block:: c++

    auto Sweep(d, u, v, q, n, a) {
      for (f : faces) {
        for (l : layers) {
          c = AdjacentCell(f, v[f] > 0 ? 0 : 1); // upwind cell
          vu(f, l) = PlicFlux(n(c, l), a(c, l), h, v, dt, d);

          cd = AdjacentCell(f, v > 0 ? 1 : 0); // downwind cell
          if (FindLayer(q, cd, q(c, l)) == kLayerNone) {
            ld = FindLayer(q, cd, kClNone);
            if (ld != kLayerNone) {
              q(cd, ld) = q(c, l);
            }
          }
        }
      }
      for (c : cells) {
        for (l : layers) {
          fm = AdjacentCell(c, 0);
          fp = AdjacentCell(c, 1);
          ds = (v[fp] - v[fm]) * dt / volume;
          vm = v(c, FindLayer(qf, fm, q(c, l)));
          vp = v(c, FindLayer(qf, fp, q(c, l)));
          dl = (vp - vm) * dt / vol;
          u(c, l) += ud[c] * ds - dl;
          if (u(c, l) == 0) {
            q(c, l) = kClNone;
          }
        }
      }
      return {u, q};
    }

Function ``Plic()`` builds the PLIC reconstruction in interfacial cells
and returns normals *n* and plane constants *a*:

.. code-block:: c++

    auto Plic() {
      for (c : cells) {
        for (l : layers) {
          uu = Stencil3(c, u, q, q(c, l));
          n(c, l) = Normal(uu);
        }
        auto navg = Avg(n, q, c);
        for (l : layers) {
          n(c, l) = (n.dot(navg) > 0 ? -navg : navg);
        }
        a(c, l) = PlaneConstant(u(c, l), n(c, l), h);
      }
      return {n, a};
    }

Normals *n* are estimated using the standard Youngs-centered scheme ``Normal()``
applied to stencil values assembled by ``Stencil3()``,
which selects values with given color (or zero if the requested color is not
found) on the *3x3x3* stencil.
After all advection sweeps,
``Recolor()`` detects new connected components
and recomputes the colors as described in :ref:`s:vofm_recolor`.

The complete implementation of the multilayer VOF
is given in :linkpath:`solver/vofm.ipp`.

.. _s:vofm_recolor:

Connected-component labeling
----------------------------

