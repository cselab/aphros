#!/usr/bin/env python

def Exp():
  # experiment Soto2018
  # Figure 2

  # SI units

  # bubble diameter
  d = 600e-6
  # bubble radius
  r = d * 0.5
  # surface tension
  sigma = 58.7e-3
  # liquid density
  rho1 = 1014.5
  # viscosity
  mu = 0.964e-3
  # capillary time (from partstr)
  T = (rho1 * r ** 3 / sigma) ** 0.5
  # Ohnesorge number
  Oh = mu / (rho1 * d * sigma) ** 0.5
  # time between snaphots (a-o)
  ts = 75e-6
  # dimensionless time between snapshots
  ts_T = ts / T

  for key in sorted(locals()):
      print("{:} \t=\t {:}".format(key, locals()[key]))
  print("\n")


Exp()
