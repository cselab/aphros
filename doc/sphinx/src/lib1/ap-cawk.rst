.. awk for CSV files

AWK for CSV files
=================

Synopsis
--------

ap.cawk [-v var=value] ['prog'] [CSV ..]


Description
-----------


Examples
--------

::

   $ cat example.csv
   x,y
   1,10
   2,20
   $ ap.cawk '{print $x + $y}' example.csv
   11
   22

::
   
   $ ap.cawk 'NR == 2 {print HEADER} {print $x + $y}' q
   x,y
   11
   22
   33
   
Source
------

.. linkpath::`deploy/scripts/cawk`
