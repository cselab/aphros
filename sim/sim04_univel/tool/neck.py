#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Read uniform grid data
# p: path
# Format:
# <nx> <ny> <nz>
# <u[0,0,0]> <u[0,0,1]> ...
# Return:
# array of shape (nx, ny, nz)
# None if file not found
def ReadPlain(p):
  if not os.path.isfile(p):
    return None
  with open(p) as f:
    ll = f.readlines()
    # shape x,y,z
    s = np.array(ll[0].split(), dtype=int)
    # shape z,y,x
    ss = tuple(reversed(s))
    # data flat
    u = np.array(ll[1].split(), dtype=float)
    # data z,y,x
    u = u.reshape(ss)
    return u

# Read scalar hdf field.
# p: path to hdf file of shape (nx, ny, nz, 1)
# Return:
# array of shape (nx, ny, nz)
# None if file not found
def ReadHdf(p):
  if not os.path.isfile(p):
    return None
  h = h5py.File(p, 'r')
  u = h['data']
  u = np.array(u)
  h.close()
  return u[:,:,:,0]

# Read file depending on extension
# .h5: hdf
# .dat: plain
def ReadExt(p):
  e = os.path.splitext(p)[1]
  if e == ".h5":
    return ReadHdf(p)
  elif e == ".dat":
    return ReadPlain(p)
  else:
    assert False


av = sys.argv

if len(av) < 4 or av[1] == '-h':
  print("""usage: ./{:} DAT X Y
Prints range of value 1 in 3D field along z at x=X y=Y.
DAT: hdf-file with extension .h5 of shape (nx,ny,nz,1)
     or plain dat file with extension .dat
X,Y: coordinates, fractions of nx,ny
Prints:
  z0 z1
z0,z1: fractions of nz, lowest face 0, highest face 1
""".format(os.path.basename(av[0])))
  exit(1)

f = av[1]
x = float(av[2])
y = float(av[3])
assert 0 <= x and x <= 1
assert 0 <= y and y <= 1

u = ReadExt(f)

nz,ny,nx = u.shape[:3]

ix = int(x * (nx - 1))
iy = int(y * (ny - 1))

uz = u[:,iy,ix]


iz = range(nz)

# index of center of uz=1
uzs = uz.sum()
izc = int((iz * uz).sum() / uzs) if uzs else 0

# height functions
hz0 = uz[:izc].sum()
hz1 = uz[izc:].sum()

z0 = (izc - hz0) / nz
z1 = (izc + hz1) / nz

print("{:} {:}".format(z0, z1))
