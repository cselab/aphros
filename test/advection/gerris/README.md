Based on 
3.2 Time-reversed VOF advection in a shear flow
from Gerris test suite

http://gerris.dalembert.upmc.fr/gerris/tests/tests/shear.html

Modifications:

* uniform mesh
* remove eps output and stat
* add vtk output


Run and convert to `.dat` suitable for `../plot.py`
  ./run
  ./convert

Plot all `u*.dat`
  ../plot.py

