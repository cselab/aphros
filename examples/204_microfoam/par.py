np = 2  # cores
nx = 80  # mesh size

from os import environ
if "FINE" in environ:
    np = 40
    nx = 160
del environ

dim = 2
conv = "exp"

bublength = 3.
Lwide = 0.5
Lnarrow = 0.25
Hwide = 0.3
Ca = 0.1

rho = 0.001

tmax = 6.
dump_field_dt = 0.05
dump_traj_dt = dump_field_dt
dumplist = ""
dumppoly = 1
