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
   $ ch.cawk '{print $x + $y}' example.csv
   11
   22
   33

Source
------

.. linkpath:: deploy/tool/cawk
