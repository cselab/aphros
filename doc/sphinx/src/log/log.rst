Surface tension as gradient
---------------------------

2018-08-30 22:08:47

**Goal**: zero integral of force

**What**: surface tension as :math:`\nabla (\kappa \alpha)` on faces,
curvature in cells copied from neighbour faces if nan.

**Result**: instability for test of single drop equilibrium,

**Data**: `log01_grad_ka <log01_grad_ka>`__

-  `grad_ka.mp4 <log01_grad_ka/grad_ka.mp4>`__: surface tension as
   :math:`\nabla (\kappa * \alpha)`
-  `k_grad_a.mp4 <log01_grad_ka/k_grad_a.mp4>`__: surface tension as
   :math:`\kappa \nabla \alpha`
-  `k_grad_a_kmean.mp4 <log01_grad_ka/k_grad_a_kmean.mp4>`__: mean
   curvature on face if both cells contain interface

Particle normal displacement
----------------------------

2018-09-02 09:41:49

**Goal**: reduce spurious flow and deformation for single drop
equilibrium

**What**: particle strings without normal displacement, position of
central particle fixed at the interface line center

**Result**: deformation of the interface greatly reduced, probably due
to stronger coupling or penalization of deformed interfaces

**Data**: `log02_dn <log02_dn>`__

-  `dn0.mp4 <log02_dn/dn0.mp4>`__: without normal displacement
-  `dn1.mp4 <log02_dn/dn1.mp4>`__: with normal displacement

``march=native``
----------------

2018-09-09 11:59:47

**Goal**: use automatic vectorization and specific optimizations

**What**: add ``-march=native`` to ``CMAKE_C_FLAGS`` and
``CMAKE_CXX_FLAGS``

**Result**: slight improvement of performance (3.60 vs 3.76 s for
``confdiff:01:assemble``) with 64 cores on Euler

**Data**: `log03_native <log03_native>`__:

-  ``out_std``: standard flags
-  ``out_native``: ``-march=native``

``max_iter``
------------

2018-09-12 14:03:27

**Goal**: choose optimal number of iterations ``max_iter``

**What**: Taylor-Green with bubbles (``sim06``) with
``max_iter``\ =1,2,3,4,5 and ``convsc``\ =\ ``cd``,\ ``quick``.

**Result**

-  qualitatively different trajectories for ``max_iter``\ =1,2;
-  convergence at ``max_iter``\ =4
-  no major difference between ``cd`` and ``quick``,
-  no difference between ``simpler``\ =0 and 1

**Data**: `log04_maxit <log04_maxit>`__

-  ``kedr.pdf``: kinetic energy dissipation rate
-  ``traj.pdf``: x-component of trajectory of bubble 2

Number of particles and strings
-------------------------------

**Goal**: check the influence of the ``part_np`` and ``part_ns`` in
curvature estimator

**What**: Free-field coalescence of bubbles based on ``sim04_univel``

**Result**

-  slight differences in shapes comparing ``np=5,0`` and ``ns=2,3``

**Data**: `log05_coalnp <log05_coalnp>`__

-  ``gen1_ns3_np5``: configuration
-  ``nx064_ns?_np?.png``: shapes at ``t=0.466`` with ``nx=64``
-  ``nx128_ns2_np9.png``: shapes at ``t=0.466`` with ``nx=128``

Number of iterations for coalescence
------------------------------------

**Goal**: check the influence of the ``min_iter`` and ``tol`` on bubble
coalescence

**What**: Free-field coalescence of bubbles based on
``sim04_univel/case/partstrcoal``

**Result**

-  negligible difference in trajectories between default and
   ``min_iter=10``

**Data**: `log06_coaliter <log06_coaliter>`__

-  ``nx*iter*``: configuration, data, movie
-  ``a.gnu``: gnuplot script

Relaxation parameters
---------------------

2019-03-21 22:48:56

**Goal**: check the influence of ``vrelax`` and ``prelax`` on bubble
coalescence

**What**: Near-wall coalescence of bubbles with
``vrelax=0.8,prelax=0.8`` and ``vrelax=0.9,prelax=1``

**Result**

-  spurious oscillations in the pressure field (``wforce``) with
   ``vrelax=0.9`` after the bubble detachment
-  no effect on x- and y-trajectories, slower rising in z after
   detachment
-  changed default parameters for ``sim04`` to ``prelax=0.8``,
   ``vrelax=0.8`` and ``min_iter=4``

**Data**: `log07_vrelax <log07_vrelax>`__

-  ``gen1``: generator of parameters for ``sim04/gen``
-  ``wforce``: force acting on ``z=0`` and ``z=1``
-  ``out``, ``stat.dat``: log and statistics

Reduced usage of ``GetCenter``
------------------------------

2019-04-27 01:19:03

**Goal**: Performance improvement.

**What**: Reduced usage of ``GetCenter`` and ``GetVectToCell``. Test of
the Taylor-Green vortex with bubbles.

**Result**: Speedup of 15%.

**Data** `log08_getcenter <log08_getcenter>`__

-  ``tg``: setup
-  ``tg0_*``: timings before
-  ``tg1_*``: timings after

Debug linear solver with ``ConvertLsCompact``
---------------------------------------------

2019-08-15 10:04:39

**Goal**: Fix ``sim12_ringgauss`` after reimplementing ``simple.ipp``
with ``ConvertLsCompact`` in ``48fa3cf0``. Iterations for pressure
correction diverged.

**What**: Bisection between ``1ffd55c5`` from ``Jan 12 2019`` and
``9cb90989`` from ``Aug 14 2019``. Complication from another bug in
``InitVort()`` which attempted to use uninitialized ``fs_`` in
``hydro.h`` (fixed by ``patch``).

**Result**: Found the problem by printing the coefficients of the linear
system before and after ``ConvertLs()`` or ``ConvertLsCompact()``. The
problem appeared from changing the order of stencil cells. In ``sim12``,
the symmetric linear solver is called twice: for initialization of
vorticity and for pressure correction. Then Hypre is apparently
initialized only once and that instance is used for both cases which
leads to the wrong order of coefficients.

**Data**: `log09_debug_linear <log09_debug_linear>`__

-  ``patch``: fix for ``InitVort()``
-  ``r``: tool to apply patch and rebuild
-  ``simple_{bad,good}.ipp``: two versions of ``simple.ipp``
-  ``min``: minimal simulation setup

Registered modules in Ubuntu
----------------------------

2020-06-03 11:34:54

Ubuntu has ``--as-needed`` by default (seen from ``gcc -dumpspecs``)
which makes the linker ignore unused libraries
and, in particular, the implementations of modules (e.g. ``init_contang.so``).

VTK merge, comparison of floats with tolerance
----------------------------------------------

2020-06-13 21:56:35

``ConvertMerge`` in ``dump/vtk.h`` mapped each ``Vect`` to a hash
and used the hash to merge closely located points together
(up to tolerance ``tol``).
This effectively split the space of ``Vect`` to partitions.
However, if two points are located close to the boundary between
two partitions, they can be arbitrarily close but have
different hash values.

The solution was to check compare against the hash
of neighboring points (from vertices of a cube)
when looking for an existing hash.

Minimal example that gave vertices of rank 1: `log10_vtkmerge <log10_vtk_merge>`__.

The opposite problem: hash collisions of distant points.
Example is in ``log10_vtkmerge/hash_collision``.
To fix, changed the map from ``hash->index`` to ``Vect->index``
(ensures exact comparison of points)
and introduced ``canonical(x)`` to get a single point from the cell of
size ``tol``.
