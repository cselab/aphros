#!/usr/bin/env python

"""
Measure the size and radius of the coalescence neck.
STDOUT:
  z0 z1 x0 x1 zc [r0 cx0 cz0 r1 cx1 cz1]
length measured in fractions of nz,
z0,z1: range of the neck
x0,x1: range in x computed at y=Y and z-center of mass
zc: z-center of mass or Z (if --z provided)
r0,r1: radius of the circular arcs
cx0,cz0,cx1,cz1: centers of the circular arcs
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
      '''create png named after DAT showing
the slice y=Y with the field and circular arcs'''
)
  p.add_argument('--th', type=float, default=0.1, help=
      '''threshold for the fitting error
      used for choosing the fitting range'''
)
  p.add_argument('--x', type=float, default=0.5, help="center of the neck is X*nz")
  p.add_argument('--y', type=float, default=0.5, help="center of the neck Y*nz")
  p.add_argument('--z', type=float, default=None, help="fixed z=Z*nz for computing x0,x1")

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
header += ['x0', 'x1', 'zc']
if fitcirc:
  header += ['r0', 'cx0', 'cz0', 'r1', 'cx1', 'cz1']

if args.header:
  print(" ".join(header))
  exit()

assert os.path.isfile(datafile), "file '%r' not found" % datafile

# output data
out = dict()

u = ReadExt(datafile)

nz,ny,nx = u.shape[:3]

# indices of neck center
ix = int(x * (nz - 1))
iy = int(y * (nz - 1))
ix = np.clip(ix, 0, nx - 1)
iy = np.clip(iy, 0, ny - 1)

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


# range in x
def GetRangeX(u, fixz=None):
  if fixz is not None:
    zc = fixz * nz
  else:
    zz = np.mgrid[0:nz,0:ny,0:nx][0]
    assert zz.ptp() == nz - 1, "invalid zz.ptp() = %r" % zz.ptp()
    # z-center of mass
    zc = (u * zz).sum() / u.sum()
  zc = np.clip(int(zc), 0, nz - 1)
  # column
  ux = u[zc,iy,:]
  # height functions
  hx0 = ux[:ix].sum()
  hx1 = ux[ix:].sum()
  # result
  x0 = (ix - hx0) / nz
  x1 = (ix + hx1) / nz

  return x0, x1, float(zc) / nz

out['x0'], out['x1'], out['zc'] = GetRangeX(u, args.z)

# slice y=Y
uxz = u[:,iy,:]


if args.fitcirc:
  import scipy.optimize as opt
  # circular arc
  def F(x, r):
    h = (np.maximum(1 - (x / r) ** 2, 0)) ** 0.5
    return (h - 1.) * r

  def GetX(nx):
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
    i1 = min(nx, ix + w + 1)
    # height from the circle
    def H(r):
      return F(xx[i0:i1] - ix, r)
    # variance
    def E(r):
      return (H(r) - hh[i0:i1]).var()
    a = opt.minimize(E, r0)
    r = a.x[0]
    return r, i0, i1, F(xx - ix, r)

  # Finds maximum range for fitting circles
  # with the error is below the threshold.
  # hh: height along x, shape (nx)
  # ix: index of center
  # th: threshold for maximum error
  # r0: initial radius
  # Returns:
  # r: optimal radius
  # i0,i1: fitting range, [ix-w,ix+w]
  # zz: height from the circle, shape (nx)
  def FindRange(hh, ix, th, r0):
    nx = hh.shape[0]
    xx = GetX(nx)
    rrprev = None
    for w in range(1,nx//2):
      rr = FitR(hh, ix, w, r0)
      r, i0, i1, zz = rr
      zz += (hh[i0:i1] - zz[i0:i1]).mean()
      e = (zz[i0:i1] -hh[i0:i1]).std()
      if e > th:
        break
      rrprev = rr
    return rrprev if rrprev is not None else rr


  # height function upper (plus)
  hhp = izc + uxz[izc:,:].sum(axis=0)
  r0 = -nx
  r, i0, i1, zz = FindRange(hhp, ix, args.th, r0)
  out['r0'] = -r / nz
  out['cx0'] = float(ix) / nz
  out['cz0'] = zz[nx // 2] / nz
  rr0 = [i0, i1, zz]

  # height function lower (minus)
  hhm = izc - uxz[:izc,:].sum(axis=0)
  r0 = nx
  r, i0, i1, zz = FindRange(hhm, ix, args.th, r0)
  out['r1'] = r / nz
  out['cx1'] = float(ix) / nz
  out['cz1'] = zz[nx // 2] / nz
  rr1 = [i0, i1, zz]

if args.plot:
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt

  plt.imshow(np.flipud(uxz), extent=(0, nx, 0, nz),
      cmap=plt.get_cmap('gray_r'), alpha=0.5)
  plt.xlim(0, nx)

  if args.fitcirc:
    xx = GetX(nx)

    i0,i1,zz = rr0
    plt.plot(xx[i0:i1], zz[i0:i1])

    i0,i1,zz = rr1
    plt.plot(xx[i0:i1], zz[i0:i1])

    plt.scatter(out['x0'] * nz, out['zc'] * nz, s=0.3)

  o = os.path.splitext(os.path.basename(datafile))[0] + ".png"
  plt.savefig(o)


# output
print(" ".join([str(out[h]) for h in header]))

