#!/bin/env python3

import math
import scipy.optimize
import sys

me = "plane.py"

def usg():
    sys.stderr.write("%s -a 'capilary constant' [-n points]\n" % me)
    sys.exit(2)

def f(r):
    return -1/math.sqrt(2) * math.acosh(math.sqrt(2)/r) + math.sqrt(2 - r**2)

a = None
t = None
n = 20
while True:
    sys.argv.pop(0)
    if len(sys.argv) and len(sys.argv[0]) > 1 and sys.argv[0][0] == '-':
        if sys.argv[0][1] == 'h':
            usg()
        elif sys.argv[0][1] == 'n':
            sys.argv.pop(0)
            n = int(sys.argv[0])
        elif sys.argv[0][1] == 'a':
            sys.argv.pop(0)
            a = float(sys.argv[0])
        elif sys.argv[0][1] == 't':
            sys.argv.pop(0)
            t = float(sys.argv[0])            
        else:
            sys.stderr.write("%s: unknown option '%s'\n" % (me, sys.argv[0]))
            sys.exit(2)
    else:
        break

if a == None:
    sys.stderr.write("%s: -a is not set\n" % me)
    sys.exit(2)

if t == None:
    sys.stderr.write("%s: -t is not set\n" % me)
    sys.exit(2)

t *= math.pi/180

L = 10
h = math.sqrt(1 - math.sin(t))
sys.stderr.write("%g\n" % h)
x0 = - f(h)

dz = h / n
for i in range(n):
    z =  dz * (i + 1/2)
    x = f(z) + x0

    x *= a
    z *= a
    print("%.16e %.16e" % (-x, z))

# L.D. Landau & E.M. Lifshitz Fluid Mechanics ( Volume 6 of A Course
# of Theoretical Physics ) Pergamon Press 1959
# p. 235
# Determine the shape of the surface of a fluid in a gravitational
# field and bounded on one side by a vertical plane wall. The angle of
# contact between the fluid and the wall is `theta' (Fig. 32).
