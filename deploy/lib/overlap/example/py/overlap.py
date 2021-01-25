from ctypes import *
import os

p = os.getenv("APHROS_PREFIX")
if p == None:
    sys.stderr.write("APHROS_PREFIX is not set")
    sys.exit(2)
p = os.path.join(p, "lib", "liboverlap.so")
l = cdll.LoadLibrary(p)

d = c_double

l.overlap_2d.restype = d
l.overlap_2d.argtypes = [d, d, d]

l.overlap_3d.restype = d
l.overlap_3d.argtypes = [d, d, d, d]

def overlap_2d(x, y, r):
    return l.overlap_2d(x, y, r)

def overlap_3d(x, y, z, r):
    return l.overlap_3d(x, y, z, r)

def Main():
    a = overlap_2d(0.1, 0.2, 0.1)
    b = overlap_3d(0.3, 0.1, 0.2, 0.1)
    print(a * 100)
    print(b * 1000 / 4 * 3)
    print(overlap_3d(10.0, 0.0, 0.0, 10))
    print(overlap_2d(10.0, 0.0, 10))

if __name__ == "__main__":
    Main()
