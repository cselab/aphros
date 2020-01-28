.. _s:embed:

Embedded boundaries
===================

The method of embedded boundaries allows
one to solve equations in complex geometries on a Cartesian mesh.
The surface of a solid body is represented by *cut cells*.
A cut cell is obtained by cutting a regular cubic cell with a plane.
Discretization of the equations with the finite volume method
remains the same, except that the discrete conservation law now includes
additional boundaries in cut cells.

The functionality related to embedded boundaries is provided in header
``embed.h``.

.. code-block:: cpp

  #include <aphros/solver/embed.h>

A new type ``FieldEmbed`` transparently combines ``FieldFace`` and
``FieldCell``.
The cell-field stores the values on embedded boundaries.
Given a mesh ``m``,
a combined field is initialized on the mesh similar to other fields
with an optional default value:

.. code-block:: cpp

  FieldEmbed<Scal> fe(m, 0); // [f]ield [e]mbed

The field can be accessed both by cell- and face-indices.
The following code traverses all internal cells and faces and assigns a
value to the field:

.. code-block:: cpp

  for (IdxCell c : m.Cells()) {
    fe[c] = 1;
  }
  for (IdxFace f : m.Faces()) {
    fe[f] = 1;
  }

A new class ``Embed<M>`` defines
new iterators and geometric routines that take into account cut cells.
The cells and faces are divided into three groups:

* regular cells (whole cells inside the domain),
* cut cells (cells crossed by the domain boundary),
* excluded cells (cells outside the domain).

The same classification applies to faces

* regular faces (whole faces inside the domain),
* cut faces (faces crossed by the domain boundary),
* excluded faces (cells outside the domain).

An instance of ``Embed<M>``, commonly named as ``eb``, is constructed from a
mesh and needs to be initialized from a level-set function :math:`\varphi`
stored in a node-field:

.. code-block:: cpp

  void F(const FieldNode<Scal>& phi) {
    auto sem = m.GetSem();
    struct {
      std::unique_ptr<Embed<M>> eb_;
    } * ctx(sem);
    if (sem()) {
      eb_.reset(new Embed<M>(m));
    }
    if (sem.Nested()) {
      eb_->Init(phi);
    }
  }

Regular cells are those for which :math:`\varphi>0` in all adjacent nodes,
excluded cells are those for which :math:`\varphi<0` in all adjacent nodes,
and other cells are cut cells.
The ``Embed<M>`` class redefines ranges ``Cells()`` and ``Faces()``
such that they include regular and cut cells but not excluded cells.
Another range ``CFaces()`` traverses only cut cells.
The name ``CFaces()`` comes from *Cell Faces*.
Together, ``Faces()`` and ``CFaces()`` include all indices
for which an embed-field can be defined.
The following code traverses non-excluded faces and cut cells
and assigns a value to the field

.. code-block:: cpp

  for (IdxCell c : eb.CFaces()) {
    fe[c] = 1;
  }
  for (IdxFace f : eb.Faces()) {
    fe[f] = 1;
  }


