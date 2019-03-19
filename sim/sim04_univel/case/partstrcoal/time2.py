#!/usr/bin/env python

def Exp():
  # experiment Thoroddsen 2005

  # SI units

  sc = 118 * 1e3   # pixels per 1 m
  # top bubble size x
  d1x = 311 / sc
  # top bubble size z
  d1z = 291 / sc
  # bottom bubble size x
  d2x = 264 / sc
  # bottom bubble size z
  d2z = 283 / sc
  b1rx = d1x * 0.5
  b1rz = d1z * 0.5
  b2rx = d2x * 0.5
  b2rz = d2z * 0.5

  # average bubble diameter
  davg = (d1x + d2x) * 0.5
  # average bubble radius
  ravg = davg * 0.5
  # surface tension
  sigma = 21.5e-3
  # liquid density
  rho1 = 789
  # viscosity
  mu = 1.2e-3
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

  # sim
  # average bubble radius
  s_ravg = 0.15
  # neck position in z
  s_neckz = 0.5
  s_b1rx = b1rx * s_ravg / ravg
  s_b1rz = b1rz * s_ravg / ravg
  s_b1cz = s_neckz + s_b1rz
  s_b2rx = b2rx * s_ravg / ravg
  s_b2rz = b2rz * s_ravg / ravg
  s_b2cz = s_neckz - s_b2rz

  for key in sorted(locals()):
      print("{:} \t=\t {:}".format(key, locals()[key]))
  print("\n")

  print('''ravg={s_ravg}
b1cz={s_b1cz}
b1rx={s_b1rx}
b1rz={s_b1rz}
b2cz={s_b2cz}
b2rx={s_b2rx}
b2rz={s_b2rz}
'''.format(**locals()))


Exp()
