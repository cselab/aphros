.. _s:bclist:

Boundary conditions
===================

Header ``func/init_u.h`` provides routines to describe boundary
conditions based on primitives.
The list of boundary conditions is a list of code blocks with primitives.

Code blocks
-----------

A code block is described by

.. includecode:: src/parse/codeblocks.h
  :struct: CodeBlock

and has the following format

.. code-block:: none

  NAME {
    CONTENT
  }

where ``NAME`` does not contain ``{`` or ``}``
and ``CONTENT`` may contain only matching ``{...}``.

Examples of code blocks

.. includecode:: test/parser/codeblocks_in
  :language: none

The following function extracts a list of blocks from a stream

.. includecode:: src/parse/codeblocks.h
  :func: ParseCodeBlocks
  :comment:

Leading and trailing whitespaces are stripped.


Groups of primitives
--------------------

Boundary conditions are selected inside groups of primitives
organized as code blocks.
The name of the block is the boundary condition to set
and the content is a list of primitives
describing a single level set function.
Affected are cut cells and faces on the domain boundary.

.. code-block:: none


  wall 0 0 0 {
    box 0 0 0 10
  }
  wall 1 0 0 {
    sphere 0.7 0.5 0.5 0.2
  }
  slipwall {
    sphere 0.3 0.5 0.5 0.2
  }
  inlet 1 0 0 {
    sphere 0.5 0.7 0.5 0.2
  }
  symm {
    sphere 0.5 0.5 0.7 0.2
    sphere 0.5 0.5 0.3 0.2
  }

The following function parses parses a file
and returns a ``MapEmbed`` object containing face and cut-cell conditions.

.. includecode:: src/func/init_bc.h
  :func: ParseGroups
  :dedent: 2
  :comment:
