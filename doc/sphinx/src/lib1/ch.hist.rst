.. histogram tool

Histogram
=========

Synopsis
--------

ch.hist [--density|--cumulative] [--bins BINS] [--range L R]

print histogram of STDIN as columns [center,density,left,right]

Examples
--------

::

   seq 10 | ch.hist --density
   seq 10 | ch.hist --cumulative

.. linkpath:: tool/hist

