#!/usr/bin/env python

from glob import glob
import re
import numpy as np
import scipy
import scipy.interpolate
import sys

def GetX(d):
    return np.loadtxt(d + "/x")

def GetY(d):
    return np.loadtxt(d + "/vy")

def GetN(d):
    m = re.findall("[a-z]*(\d*)", d)[0]
    return int(m)

dd = glob("ch*")

# dict {n:d}
ddn = dict()
for d in dd:
    ddn[GetN(d)] = d

# finest mesh
d = ddn[max(ddn)]
nm = GetN(d)
xm = GetX(d)
print("finest nx={:}".format(nm))

# dict {n:y}
yy = dict()

# interpolate to finest
for n in sorted(ddn):
    d = ddn[n]
    print(d)
    x = GetX(d)
    y = GetY(d)
    f = scipy.interpolate.interp1d(x, y,
            kind='nearest', fill_value='extrapolate')
    yy[n] = f(xm)

# write error
with open('er', 'w') as o:
    o.write("nx e1 e2 em\n")
    for n in sorted(ddn):
        ym = yy[nm]
        y = yy[n]
        e1 = (abs(y - ym)).mean()
        e2 = ((y - ym) ** 2).mean() ** 0.5
        em = abs(y - ym).max()

        o.write("{:} {:} {:} {:}\n".format(n, e1, e2, em))

# write difference profile
for n in sorted(ddn)[:-1]:
    d = ddn[n]
    ym = yy[nm]
    y = yy[n]
    np.savetxt(d + "/dvy", np.array((xm, y - ym)).T)
