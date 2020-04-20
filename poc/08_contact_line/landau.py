#!/bin/env python3

import math
import scipy.optimize
import sys

me = "lambda.py"

def f(r, z, c):
    return c * math.cosh(z/c)

n = 20
d = 0.5
r = 1.0
while True:
    sys.argv.pop(0)
    if len(sys.argv) and len(sys.argv[0]) > 1 and sys.argv[0][0] == '-':
        if sys.argv[0][1] == 'n':
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

eps = 1e-2
c0 = eps
g = lambda c : f(r, d, c) - r
res = scipy.optimize.minimize(g, c0, bounds = ((eps, None),))
if not res.success:
    sys.stderr.write("%s: minimize failed\n" % me)
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
dz = d / (n - 1)
pos = [ ]
neg = [ ]
for i in range(n):
    z = dz * i
    r = g(z)
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
