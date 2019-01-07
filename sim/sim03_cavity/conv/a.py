#!/usr/bin/env python

from glob import glob
import re
import numpy as np
import scipy
import scipy.interpolate

dd = sorted(glob("ch*"))

def GetX(d):
    return np.loadtxt(d + "/x")

def GetY(d):
    return np.loadtxt(d + "/vy")

def GetN(d):
    return int(re.findall("ch(\d*)", d)[0])

# finest mesh
d = dd[-1]
nm = GetN(d)
xm = GetX(d)

# (nx:y)
yy = dict()

for d in dd:
    print(d)
    n = GetN(d)
    x = GetX(d)
    y = GetY(d)
    f = scipy.interpolate.interp1d(x, y,
            kind='nearest', fill_value='extrapolate')
    yy[n] = f(xm)

# write error
with open('er', 'w') as o:
    o.write("nx e1 e2 em\n")
    for n in yy:
        ym = yy[nm]
        y = yy[n]
        e1 = (abs(y - ym)).mean()
        e2 = ((y - ym) ** 2).mean() ** 0.5
        em = abs(y - ym).max()

        o.write("{:} {:} {:} {:}\n".format(n, e1, e2, em))

# write difference profile
for d in dd:
    n = GetN(d)
    ym = yy[nm]
    y = yy[n]
    np.savetxt(d + "/dvy", np.array((xm, y - ym)).T)
