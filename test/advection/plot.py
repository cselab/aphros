#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# p -- path
def Read(p):
    with open(p) as f:
        ll = f.readlines()
        # size
        s = np.array(ll[0].split(), dtype=int)
        # data
        u = np.array(ll[1].split(), dtype=float)
        return u.reshape(s)

# u -- numpy array (2d or 3d with shape[2]==1)
def Get2d(u):
    if len(u.shape) == 2:
        return u
    else:
        s = u.shape
        assert len(s) == 3 and s[2] == 1
        return u.reshape((s[0], s[1]))

# u -- 2d numpy array
# fn -- output file name
def Plot(u, fn):
    plt.imshow(np.flipud(u), extent=(0, 1, 0, 1), interpolation='nearest')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    print(fn)
    plt.savefig(fn, dpi=300)
    plt.close()

for i in range(3):
    b = "u" + str(i)
    u = Read(b + ".dat")
    u = Get2d(u)
    Plot(u, b + ".pdf")
