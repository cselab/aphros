import numpy
import os

# Read uniform grid data
# p: path
# Format:
# <nx> <ny> <nz>
# <u[0,0,0]> <u[0,0,1]> ...
# Return:
# array of shape (nx, ny, nz)
# None if file not found
def ReadPlain(fn):
  if not os.path.isfile(fn):
    return None
  with open(fn) as f:
    ll = f.readlines()
    # shape x,y,z
    s = numpy.array(ll[0].split(), dtype=int)
    # shape z,y,x
    ss = tuple(reversed(s))
    # data flat
    u = numpy.array(ll[1].split(), dtype=float)
    # data z,y,x
    u = u.reshape(ss)
    return u
