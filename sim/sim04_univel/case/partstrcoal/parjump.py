#!/usr/bin/env python

def Exp():
  # experiment Thoroddsen 2005

  # SI units

  # equivalent bubble radius
  ravg = 0.79e-3
  # equivalent bubble diameter
  davg = 2 * ravg
  # surface tension
  sigma = 72e-3
  # liquid density
  rho1 = 1000
  # viscosity
  mu = 1e-3
  # capillary time
  T = (rho1 * ravg ** 3 / sigma) ** 0.5
  # Ohnesorge number
  Oh = mu / (rho1 * davg * sigma) ** 0.5
  # time between snaphots in Fig. 2 (top row)
  ts = 50e-6
  # max time
  tm = 2200e-6
  # dimensionless time between snapshots
  ts_T = ts / T
  # number of snapshots
  tm_ts = tm / ts
  # dimensionless max time
  tm_T = tm / T
  # gravity
  g = 9.8
  # Eotvos number
  Eo = rho1 * g * davg ** 2 / sigma

  for key in sorted(locals()):
      print("{:} \t=\t {:}".format(key, locals()[key]))
  print("\n")

Exp()
