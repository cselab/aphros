.. awk for CSV files

AWK for CSV files
=================

Synopsis
--------

ch.cawk [-v var=value] ['prog'] [file ..]


Examples
--------

::

   $ cat example.csv
   x,y
   1,10
   2,20
   3,30
   4,40

   $ ch.cawk '{print $x + $y}' example.csv
   11
   22
   33
   44

.. linkpath:: deploy/tool/cawk
