# Cubism-hydro

Distributed version of **hydro** based on **Cubism**.

## Clone

    git clone git@gitlab.ethz.ch:kpetr/Cubism-hydro.git

## Build

### Hypre

    cd lib/build
    ./all

### Cubism-hydro

    cd src
    ./conf
    make

    (on daint use `build.daint` instead)

## Run
  
    cd sim/sim01
    ./job

## Other repos

### hydro

    git clone https://github.com/divfree/hydro.git 

### mfer

    git clone git@gitlab.ethz.ch:mavt-cse/mfer.git

## Intro

### Solver

Finite volume solver for multiphase flows.

Fluid solver based on SIMPLE 
(Semi-Implicit Method for Pressure Linked Equations) 
method coupled
with TVD Superbee advection scheme for concentration (diffuse interface method).

Advantages compared to Gerris:

* moving mesh: ability to solve in the frame of reference of bubble 
which reduces numerical diffusion 

* conservation of mass up to machine precision 
(depending of the accuracy of linear solvers)

* potentially high performance of CUBISM framework
  - cache-aware (computation by blocks)
  - hiding of communication time (inner blocks processed first)
  - distributed MPI I/O 

* extendability to multiphysics (easy to add more fluid components or phenomena)


Design following CUBISM framework 

* layers:   
  - core: computation in a single block
  - node: processing of blocks on a single node 
  - cluster: communication between nodes

* ability to run on a single node without MPI

### Simulations:

* hydrodynamic interaction of bubbles
* bubble coalescence
* bubble growth of surfaces


  
