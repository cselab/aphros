# Cubism-hydro

Distributed version of **hydro** based on **Cubism**.

## Clone

    git clone git@gitlab.ethz.ch:mavt-cse/Cubism-hydro.git
    cd Cubism-hydro
    git config pull.rebase true

## Build and install

*   Follow instructions from `deploy/README.md` to
prepare environment and install dependencies.
*   Configure, build and install

     ```
     cd src
     ./conf
     # ./conf -DOVERLAP=0  # to build without overlap and Eigen
     # ./conf -CMAKE_SKIP_INSTALL_ALL_DEPENDENCY=1 # to prevent build on
     install
     make
     # make target=mfer    # to only build ch.mfer
     ```

*   Run tests

     ```
     make test
     ```

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
(depending on the accuracy of linear solvers)

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
* bubble growth on surfaces


  
