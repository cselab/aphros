#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

av = sys.argv

if len(av) < 4 or av[1] == '-h':
  print("""usage: ./{:} HDF0 X Y
Prints range of value 1 in field HDF along z at x=X y=Y.
HDF: hdf-file of shape (nx,ny,nz,1)
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

h = h5py.File(f, 'r')
u = h['data']
u = np.array(u)
h.close()

nz,ny,nx = u.shape[:3]


ix = int(x * (nx - 1))
iy = int(y * (ny - 1))

uz = u[:,iy,ix,0]


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
