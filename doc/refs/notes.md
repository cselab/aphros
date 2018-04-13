# Popinet2009

Popinet, 2009: An accurate adaptive solver for 
surface-tension-driven interfacial flows.

## Geometrical VOF scheme

Two steps:

* interface reconstruction
* geometrical flux computation and interface advection

### Interface reconstruction

The interface is represented in each cell by a line.
Given normal vector, location of the interface and volume fraction are 
in correspondence.

Conversion is done with routines by [52] Scardovelli2000.

# Reconstruction

## Volume fraction and interface location:

[52] Scardovelli, Zaleski, 2000: Analytical relations connecting linear 
interfaces and volume fractions in rectangular grids.

## Normal estimation

[4] Scardovelli, Zaleski, 1999: Direct numerical simulation of free-surface 
and interfacial flow.

[53] Scardovelli, Zaleski, 2003: Interface reconstruction with least-square
fit and split eulerian-lagrangian advection

[54] Aulisa, Manservisi, Scardovelli, Zaleski, 2007:
Interface reconstruction with least-squares fit and split
advection in three-dimensional Cartesian geometry.

# Aulisa2007

Aulisa, Manservisi, Scardovelli, Zaleski, 2007:
Interface reconstruction with least-squares fit and split
advection in three-dimensional Cartesian geometry.

## Intro

_Consistency_ means $0 \leq C \leq 1$. 

