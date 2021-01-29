.. _s:comm:

Communication
=============

The grid is decomposed into subdomains, one for each MPI rank.
Each subdomain is further divided into cubic blocks of fixed size.
The program is built following the SPMD (single program multiple data) model.
Each MPI rank processes its blocks sequentially.
The code is divided into *stages*
which are executed in parallel on all blocks
with synchronization after every stage.
This implements *coroutines*.
The code can issue communication requests
which are collected during each stage
and executed before transferring to the next stage.

The following function consists of two stages:

.. code-block:: cpp

  void F(Mesh& m) {
    auto sem = m.GetSem();
    struct {
      FieldCell<Scal> fc;
    }* ctx(sem);

    if (sem()) {
      fc.Reinit(m, m.GetId());
      m.Comm(&fc);
    }

    if (sem()) {
      fc.Reinit(m, 0);
    }
  }

Mesh `m` provides functionality for parallel communication.
Pointer `ctx` will be initialized to point at an instance
of the user-defined `struct` persistent across stages 
of one call of the function.
The instance is created before the execution of the first stage
and destroyed after the last stage.

Communication requests include:

.. includecode:: src/geom/mesh.h
  :func: Comm

.. includecode:: src/geom/mesh.h
  :func: Reduce

.. includecode:: src/geom/mesh.h
  :func: GatherToLead

.. includecode:: src/geom/mesh.h
  :func: Dump
