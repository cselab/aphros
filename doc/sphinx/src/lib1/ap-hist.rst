.. histogram tool

Histogram
=========

Synopsis
--------

ap.hist [--density|--cumulative] [--bins BINS] [--range L R]

print histogram of STDIN as columns [center,density,left,right]

Examples
--------

::

   seq 10 | ap.hist --density
   seq 10 | ap.hist --cumulative


Source
------

:linkpath:`scripts/hist`
