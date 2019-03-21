## Surface tension as gradient

2018-08-30 22:08:47

**Goal**: 
zero integral of force

**What**: 
surface tension as $`\nabla (\kappa \alpha)`$ on faces,
curvature in cells copied from neighbour faces if nan.

**Result**: 
instability for test of single drop equilibrium, 

**Data**:
[`log01_grad_ka`](log01_grad_ka)

- [`grad_ka.mp4`](log01_grad_ka/grad_ka.mp4): 
surface tension as $`\nabla (\kappa * \alpha)`$
- [`k_grad_a.mp4`](log01_grad_ka/k_grad_a.mp4):
surface tension as $`\kappa \nabla \alpha`$
- [`k_grad_a_kmean.mp4`](log01_grad_ka/k_grad_a_kmean.mp4): 
mean curvature on face if both cells contain interface


## Particle normal displacement 

2018-09-02 09:41:49

**Goal**:
reduce spurious flow and deformation for single drop equilibrium

**What**:
particle strings without normal displacement,
position of central particle fixed at the interface line center

**Result**:
deformation of the interface greatly reduced,
probably due to stronger coupling or penalization of deformed interfaces

**Data**:
[`log02_dn`](log02_dn)

- [`dn0.mp4`](log02_dn/dn0.mp4): without normal displacement
- [`dn1.mp4`](log02_dn/dn1.mp4): with normal displacement


## `march=native`

2018-09-09 11:59:47

**Goal**:
use automatic vectorization and specific optimizations

**What**:
add `-march=native` to `CMAKE_C_FLAGS` and `CMAKE_CXX_FLAGS`

**Result**:
slight improvement of performance (3.60 vs 3.76 s for `confdiff:01:assemble`)
with 64 cores on Euler

**Data**: 
[`log03_native`](log03_native): 

- `out_std`: standard flags
- `out_native`: `-march=native`

## `max_iter`

2018-09-12 14:03:27

**Goal**:
choose optimal number of iterations `max_iter`

**What**:
Taylor-Green with bubbles (`sim06`) with `max_iter`=1,2,3,4,5
and `convsc`=`cd`,`quick`.

**Result**

- qualitatively different trajectories for `max_iter`=1,2;
- convergence at `max_iter`=4
- no major difference between `cd` and `quick`,
- no difference between `simpler`=0 and 1

**Data**:
[`log04_maxit`](log04_maxit)

- `kedr.pdf`: kinetic energy dissipation rate
- `traj.pdf`: x-component of trajectory of bubble 2

## Number of particles and strings

**Goal**:
check the influence of the `part_np` and `part_ns` in curvature estimator

**What**:
Free-field coalescence of bubbles based on `sim04_univel`

**Result**

- slight differences in shapes comparing `np=5,0` and `ns=2,3`

**Data**:
[`log05_coalnp`](log05_coalnp)

- `gen1_ns3_np5`: configuration
- `nx064_ns?_np?.png`: shapes at `t=0.466` with `nx=64`
- `nx128_ns2_np9.png`: shapes at `t=0.466` with `nx=128`

## Number of iterations for coalescence

**Goal**:
check the influence of the `min_iter` and `tol` on bubble coalescence

**What**:
Free-field coalescence of bubbles based on `sim04_univel/case/partstrcoal`

**Result**

- negligible difference in trajectories between default and `min_iter=10`

**Data**:
[`log06_coaliter`](log06_coaliter)

- `nx*iter*`: configuration, data, movie
- `a.gnu`: gnuplot script

## Relaxation parameters 

2019-03-21 22:48:56

**Goal**:
check the influence of `vrelax` and `prelax` on bubble coalescence

**What**:
Near-wall coalescence of bubbles with `vrelax=0.8,prelax=0.8` 
and `vrelax=0.9,prelax=1`

**Result**

- spurious oscillations in the pressure field (`wforce`) with `vrelax=0.9`
  after the bubble detachment
- no effect on x- and y-trajectories, slower rising in z after detachment
- changed default parameters for `sim04` to 
  `prelax=0.8`, `vrelax=0.8` and `min_iter=4`

**Data**:
[`log07_vrelax`](log07_vrelax)

- `gen1`: generator of parameters for `sim04/gen`
- `wforce`: force acting on `z=0` and `z=1` 
- `out`, `stat.dat`: log and statistics
