Introduction
============

**Aphrós** is a finite volume solver for incompressible multiphase flows with
surface tension. Key features:

- implementation in C++11
- abstractions for mesh elements (cells, faces and nodes),
  range-based loops over them
  and the corresponding data fields (cell-, face- and node-based fields)
- convenient and fast development (no need to write x,y,z loops)
  using only standard features of the language
  (without code generators or domain-specific languages)
- scalability to thousands of compute nodes using MPI/OpenMP
  thanks to the Cubism library for distributed computing on structured grids
- coroutines to enable encapsulation in the block-wise processing framework
- individual solvers can be used separately as regular classes or functions
- no Python bindings
- fluid solver based on SIMPLE or Chorin's projections
- conservative split advection solver based on PLIC
- method of particles for curvature estimation that outperforms
  standard techniques at low resolutions
  which is well-suited for simulation with many small bubbles and drops
- multilayer volume-of-fluid scheme for coalescence prevention
  which has the computational cost that does not depend on the number
  of bubbles in the simulation and therefore
  can describe processes such as generation of foam
