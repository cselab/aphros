.. _s:bclist:

Boundary conditions
===================

Header ``func/bclist.h`` provides routines to describe boundary
conditions based on primitives.

The list of boundary conditions is a list of code blocks with primitives.
A code block is described by

.. includecode:: src/func/primlist.h
  :struct: CodeBlock

and has the following format

::
  NAME {
    CONTENT
  }

where ``NAME`` does not contain ``{`` or ``}``
and ``CONTENT`` may contain only matching ``{...}``.

Examples of code blocks

::
  name {
    content
  }

  name with spaces {
    content { inside braces }
  }

  name { content }

The following function extracts a list of blocks from a stream

.. includecode:: src/func/primlist.h
  :func: ParseCodeBlocks
  :comment:

Leading and trailing whitespaces are stripped.
