#!/usr/bin/env python

"""
Measure the size and curvature of the coalescence neck.
STDOUT:
  z0 z1 [k0 x0 t1 k1 x1 y1]
where (length measured in fractions of nz)
z0,z1: range of the neck
k0,k1: curvature of the circular arcs
x0,y0,x1,y1: centers of the circular arcs
"""

import numpy as np
import sys
import os
import argparse

def Parse():
  p = argparse.ArgumentParser(description=__doc__,
      formatter_class=argparse.RawTextHelpFormatter)
  p.add_argument('--header',  action='store_true', help=
      '''print the data header and exit'''
)
  p.add_argument('datafile', nargs="?", type=str, help=
      '''path to hdf5 with extention .h5 of shape (nz,ny,nx,1)
or plain data with extention .dat of shape (nz,ny,nx);
pass 'header' to only print the output header''')
  p.add_argument('--fitcirc',  action='store_true', help=
      '''fit circular arcs to the neck'''
)
  p.add_argument('--plot', action='store_true', help=
      '''create pdf named after DAT showing
the slice y=Y with the field and circular arcs'''
)
  p.add_argument('--th', default=0.1, help=
      '''threshold for the fitting error
      used for choosing the fitting range'''
)
  p.add_argument('--x', default=0.5, help="center of the neck is X*nz")
  p.add_argument('--y', default=0.5, help="center of the neck Y*nz")

  return p.parse_args()

args = Parse()

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
  import h5py
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


datafile = args.datafile
x = args.x
y = args.y
fitcirc = args.fitcirc

header = ['z0', 'z1']
if fitcirc:
  header += ['k0', 'x0', 'y0', 'k1', 'x1', 'y1']

if args.header:
  print(" ".join(header))
  exit()

# output data
out = dict()

u = ReadExt(datafile)

nz,ny,nx = u.shape[:3]

assert 0 <= x and x <= nx / nz
assert 0 <= y and y <= ny / nz

# indices of neck center
ix = int(x * (nz - 1))
iy = int(y * (nz - 1))

# column through neck center
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

out['z0'] = z0
out['z1'] = z1


if args.fitcirc:
  import scipy.optimize as opt

  # slice y=Y
  uxz = u[:,iy,:]

  # circular arc
  def F(x, r):
    h = (np.maximum(1 - (x / r) ** 2, 0)) ** 0.5
    return (h - 1.) * r

  def GetX():
    return np.arange(nx)

  # Fits circle to height by minimizing the variance.
  # hh: height along x, shape (nx)
  # ix: index of center
  # w: range half-width, range is [ix-w,ix+w]
  # r0: initial radius
  # Returns:
  # r: optimal radius
  # i0,i1: fitting range
  # zz: height from the circle, shape (nx)
  def FitR(hh, ix, w, r0):
    nx = hh.shape[0]
    xx = GetX(nx)
    i0 = max(0, ix - w)
    i1 = min(n, ix + w + 1)
    # height from the circle
    def F(r):
      zz = F(xx[i0:i1] - ix, r)
    # variance
    def E(r):
      return (F(r) - hh[i0:i1]).var()
    a = opt.minimize(E, r0)
    r = a.x
    return r, i0, i1, F(xx - ix, r)

  # Finds maximum range for fitting circles
  # with the error is below the threshold.
  # hh: height along x, shape (nx)
  # ix: index of center
  # th: threshold for maximum error
  # Returns:
  # r: optimal radius
  # w: fitting range half-width
  # i0,i1: fitting range
  # zz: height from the circle, shape (nx)
  def FindRange(hh, th):
    nx = hh.shape[0]
    xx = GetX(nx)
    for w in range(1,nx//2):
      r, i0, i1, zz = FitR(hh, ix, w, -nx)
      zz += (hh[i0:i1] - zz[i0:i1]).mean()
      e = (zz[i0:i1] -hh[i0:i1]).std()
      print("w={:} e={:}".format(w, e))
      if e > th:
        break
      return

  out['k0'] =
  # height function upper (plus)
  hhp = izc + uxz[izc:,:].sum(axis=0)

  # height function lower (minus)
  hhm = izc - uxz[:izc,:].sum(axis=0)


  # threshold for std
  th = 0.1
  for w in range(1,30):
    r, i0, i1 = FitR(hhm, ix, w, nx)
    zz = F(xx - ix, r)
    zz += (hhm[i0:i1] - zz[i0:i1]).mean()
    e = (zz[i0:i1] - hhm[i0:i1]).std()
    print("w={:} e={:}".format(w, e))
    if e > th:
      break



o = os.path.splitext(os.path.basename(f))[0] + ".pdf"
plt.savefig(o)

print("{:} {:}".format(z0, z1))

if args.plot:
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt

  plt.imshow(np.flipud(uxz), extent=(0, nx, 0, nz),
      cmap=plt.get_cmap('Blues'))
  plt.xlim(0, nx)

  plt.plot(xx[i0:i1], zz[i0:i1])

  plt.plot(xx[i0:i1], zz[i0:i1])
