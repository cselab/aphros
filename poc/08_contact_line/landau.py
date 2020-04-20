#!/bin/env python3

import math
import scipy.optimize
import sys

me = "lambda.py"

def usg():
    sys.stderr.write("%s -d distance -r radious [-n points]\n" % me)
    sys.exit(2)

def f(r, z, c):
    return c * math.cosh(z/c)

def df(r, z, c):
    r = z / c
    return math.cosh(r) - r * math.sinh(r)

n = 20
d = None
r = None
while True:
    sys.argv.pop(0)
    if len(sys.argv) and len(sys.argv[0]) > 1 and sys.argv[0][0] == '-':
        if sys.argv[0][1] == 'h':
            usg()
        elif sys.argv[0][1] == 'n':
            sys.argv.pop(0)
            n = int(sys.argv[0])
        elif sys.argv[0][1] == 'd':
            sys.argv.pop(0)
            d = float(sys.argv[0])
        elif sys.argv[0][1] == 'r':
            sys.argv.pop(0)
            r = float(sys.argv[0])
        else:
            sys.stderr.write("%s: unknown option '%s'\n" % (me, sys.argv[0]))
            sys.exit(2)
    else:
        break

if d == None:
    sys.stderr.write("%s: -d is not set\n" % me)
    sys.exit(2)

if r == None:
    sys.stderr.write("%s: -r is not set\n" % me)
    sys.exit(2)

r0 = r
r /= r0
d /= r0

eps = 1e-2
c0 = 2 * eps
g = lambda c : f(r, d / 2, c) - r
dg = lambda c : df(r, d / 2, c)
res = scipy.optimize.minimize(g, c0, jac = dg, bounds = ((eps, None),), method = 'L-BFGS-B')
if not res.success:
    sys.stderr.write("%s: minimize failed (%s)\n" % (me, res.status))
    sys.exit(2)
cmin = res.x
if res.fun > 0:
    sys.stderr.write("%s: there is no root (cmin = %g, res.fun = %g)\n" % (me, cmin, res.fun))
    sys.exit(2)

res = scipy.optimize.root_scalar(g, method = "bisect", bracket = [eps, cmin])
if not res.converged:
    sys.stderr.write("%s: root_scalar failed\n")
    sys.exit(2)
c0 = res.root
g = lambda z : f(r, z, c0)
dz = d / 2 / (n - 1)
pos = [ ]
neg = [ ]
for i in range(n):
    z = r0 * dz * i
    r = r0 * g(z)
    pos.append((z, r))
    pos.append((-z, r))
    neg.append((z, -r))
    neg.append((-z, -r))
pos.sort()
neg.sort()
for z, r in pos:
    print("%.16e %.16e" % (r, z))
print()
for z, r in neg:
    print("%.16e %.16e" % (r, z))

# L.D. Landau & E.M. Lifshitz Fluid Mechanics ( Volume 6 of A Course
# of Theoretical Physics ) Pergamon Press 1959
# p. 234
# Problem 1. Determine the shape of a film of liquid supported on two circular frames
# with their centres on a line perpendicular to their planes, which are parallel; Fig. 31 shows a
# cross-section of the film.
