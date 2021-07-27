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

.. code-block:: cpp

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

.. code-block:: cpp

    auto Sweep(d, u, v, q, n, a) {
      // v: mixture flux
      // w: phase flux
      // qf: phase color from upwind cell
      for (f : faces) {
        for (l : layers) {
          c = AdjacentCell(f, v[f] > 0 ? 0 : 1); // upwind cell
          w(f, l) = PlicFlux(n(c, l), a(c, l), h, v, dt, d);
          qf(f, l) = q(c, l);

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
          wm = w(c, FindLayer(qf, fm, q(c, l)));
          wp = w(c, FindLayer(qf, fp, q(c, l)));
          dl = (wp - wm) * dt / vol;
          u(c, l) += udiv(c, l) * ds - dl;
          if (u(c, l) == 0) {
            q(c, l) = kClNone;
          }
        }
      }
      return {u, q};
    }

Function ``Plic()`` builds the PLIC reconstruction in interfacial cells
and returns normals *n* and plane constants *a*:

.. code-block:: cpp

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

After each advection step, the colors need to be updated to detect new
connected components.  The following function initializes the color field with
unique values, then iteratively joins the colors from adjacent cells by taking
the minimal color until equilibration.

.. code-block:: cpp

  auto Recolor(u, q) {
    q_new = InitUnique(u, q);

    changed = true;
    while (changed) {
      changed = false;
      q_new = RecolorCorners(u, q, q_new);
      for (c : cells) {
        for (l : layers) {
          for (cn : stencil3(c)) {
            ln = FindLayer(q, cn, q(c, l));
            if (ln != kLayerNone) {
              if (q_new(cn, ln) < q_new(c, l)) {
                changed = true;
                q_new(c, l) = q_new(cn, ln);
              }
            }
          }
        }
      }
    }
    return q_new;
  }

To reduce the number of iterations, function ``RecolorCorners()``
runs the same algorithm for corners from cubic subdomains.
The colors propagate through the domain faster.


.. code-block:: cpp

  auto RecolorCorners(u, q, q_new) {
    // map: mapping from old color to new color
    for (c : corners of subdomains) {
      for (l : layers) {
        for (size_t d : {0, 1, 2}) {
          cn = AdjacentCell(c, d);
          ln = FindLayer(q, cm, q(c, l));
          q1 = q_new(c, l);
          q2 = q_new(cn, ln);
          map[max(q1, q2)] = min(q1, q2);
        }
      }
    }

    changed = true;
    while (changed) {
      changed = false;
      for (q1 in map) {
        if (map[q1] in map) {
          map[q1] = map[map[q1]];
          changed = true;
        }
      }
    }

    for (f : faces) {
      cm = AdjacentCell(f, 0);
      cp = AdjacentCell(f, 1);
      for (l : layers) {
        q_new(cm, l) = map[q(cm, l)];
        q_new(cp, l) = map[q(cp, l)];
      }
    }
    return q_new;
  }
