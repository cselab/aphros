#!/usr/bin/env python

"""
Print slice in plain format to STDOUT.
--x,y,z specify the location of the slice origin.
If None, select all indices.
length measured in factions of nz.
"""

import numpy as np
import sys
import os
import argparse

def Parse():
  p = argparse.ArgumentParser(description=__doc__,
      formatter_class=argparse.RawTextHelpFormatter)
  p.add_argument('datafile', type=str, help=
      '''path to hdf5 with extention .h5 of shape (nz,ny,nx,1)
or plain data with extention .dat of shape (nz,ny,nx)''')
  p.add_argument('-x', default="", help="origin x")
  p.add_argument('-y', default="", help="origin y")
  p.add_argument('-z', default="", help="origin z")

  return p.parse_args()

args = Parse()

# Read uniform grid data
# p: path
# Format:
# <nx> <ny> <nz>
# <u[z,y,x]> ...
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

# Writes array in plain format:
# <nx> <ny> <nz>
# <u[z,x,y]> ...
def WritePlain(f, u):
  h = " ".join([str(n) for n in reversed(u.shape)])
  if not h: h = "1"
  f.write(h + '\n')
  np.savetxt(f, u.flatten(), newline='', fmt='%.16g ')


datafile = args.datafile
x = args.x
y = args.y
z = args.z

x = None if x == "" else float(x)
y = None if y == "" else float(y)
z = None if z == "" else float(z)

assert os.path.isfile(datafile), "file '%r' not found" % datafile

u = ReadExt(datafile)

nz,ny,nx = u.shape[:3]

# indices of origin
ix = np.clip(int(x * (nz - 1)), 0, nx - 1) if x is not None else None
iy = np.clip(int(y * (nz - 1)), 0, ny - 1) if y is not None else None
iz = np.clip(int(z * (nz - 1)), 0, nz - 1) if z is not None else None

# slice
sx = ix if ix is not None else ":"
sy = iy if iy is not None else ":"
sz = iz if iz is not None else ":"
c = "u[{:},{:},{:}]".format(sx, sy, sz)
uo = eval(c)

WritePlain(sys.stdout, uo)
