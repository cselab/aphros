#!/usr/bin/env python


import numpy as np
import sys

av = sys.argv
if len(av) != 1:
    print('''usage: {:}
STDIN: 2D array u in plain format
STDOUT: coordinates of points for which 0 < u < 1,
        fractions of ny, reflected along x and shifted to zero mean
'''.format(av[0]))
    exit(1)

f = sys.stdin

h = f.readline()

nx,ny = map(int ,h.split())
u = np.loadtxt(f)
u = u.reshape((ny, nx))

ii = np.where((u > 0) & (u < 1))

x = ii[1] / ny
y = ii[0] / ny

x = np.hstack((x, 1. - x))
y = np.hstack((y, y))

x -= x.mean()
y -= y.mean()

o = sys.stdout
np.savetxt(o, np.array((x, y)).T)

