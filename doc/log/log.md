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



