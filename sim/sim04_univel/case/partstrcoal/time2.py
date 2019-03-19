#!/usr/bin/env python

def Exp():
  # experiment Thoroddsen 2005

  # SI units

  # top bubble radius
  rt = 1.59e-3
  # bottom bubble radius
  rb = 1.01e-3
  # bottom to top ratio
  rbt = rb / rt
  # top bubble diameter
  dt = rt * 2
  # surface tension
  sigma = 21.5e-3
  # liquid density
  rho1 = 789
  # viscosity
  mu = 1.2e-3
  # capillary time
  T = (rho1 * rt ** 3 / sigma) ** 0.5
  # Ohnesorge number
  Oh = mu / (rho1 * dt * sigma) ** 0.5
  # time between snaphots in Fig. 2 (top row)
  ts = 50e-6
  # dimensionless time between snapshots
  ts_T = ts / T

  sc = 0.001
  # top bubble size x
  d1x = 311 * sc
  # top bubble size z
  d1z = 291 * sc
  d2x = 264 * sc
  d2z = 283 * sc
  neckz = 0.5
  b1cz = neckz + d1z * 0.5
  b1rx = d1x * 0.5
  b1rz = d1z * 0.5
  b2cz = neckz - d2z * 0.5
  b2rx = d2x * 0.5
  b2rz = d2z * 0.5

  for key in sorted(locals()):
      print("{:} \t=\t {:}".format(key, locals()[key]))
  print("\n")


Exp()
