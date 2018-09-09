2018-08-30 22:08:47

# Surface tension as gradient

Goal: 
  zero integral of force

What: 
  surface tension as $$\nabla (\kappa \alpha)$$ on faces,
  curvature in cells copied from neighbour faces if nan.

Result: 
  instability for test of single drop equilibrium, 

Data:
    log01_grad_ka:

    grad_ka.mp4: surface tension as $$\nabla (\kappa * \alpha)$$
    k_grad_a.mp4: surface tension as $$\kappa \nabla \alpha$$
    k_grad_a_kmean.mp4: mean curvature on face if both cells contain interface



2018-09-02 09:41:49

# Normal displacement from center

Goal:
  reduce spurious flow and deformation for single drop equilibrium

What:
  particle strings without normal displacement,
  position of central particle fixed at the interface line center

Result:
  deformation of the interface greatly reduced,
  probably due to stronger coupling or penalization of deformed interfaces

Data:
    log02_dn:

    dn0.mp4: without normal displacement
    dn1.mp4: with normal displacement



2018-09-09 11:59:47

# march=native

Goal:
  use automatic vectorization and specific optimizations

What:
  add `-march=native` to CMAKE_C_FLAGS and CMAKE_CXX_FLAGS

Result:
  slight improvement of performance (3.76 vs 3.60 s for confdiff:01:assemble),
  64 cores on Euler.

Data:
    log03_native: 
    out_std: standard flags
    out_native: -march=native
