from ctypes import *
import os

p = os.getenv("CHPREFIX")
if p == None:
    sys.stderr.write("CHPREFIX is not set")
    sys.exit(2)
p = os.path.join(p, "lib", "libvof.so")

l = cdll.LoadLibrary(p)

d = c_double
l.vof_cylinder.restype = d
l.vof_cylinder.argtypes = [d] * 7

def cylinder(x, y, z, r, u, v, w):
    return l.vof_cylinder(x, y, z, r, u, v, w)
