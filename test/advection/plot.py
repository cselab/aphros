#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os

# Read uniform grid data
# p -- path
# Format:
# <Nx> <Ny> <Nz>
# <u[0,0,0]> <u[1,0,0]> ...
# Return:
# array of shape (Nx, Ny, Nz)
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
# po -- output path
def Plot(u, po):
    u = np.clip(u, 0., 1.)
    fig, ax = plt.subplots(figsize=(5,5))
    ax.imshow(np.flipud(u),
               extent=(0, 1, 0, 1),
               interpolation='nearest')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(po, dpi=300)
    plt.close()



pp = sorted(glob.glob("u*.dat"))
for p in pp:
    u = Read(p)
    u = Get2d(u)
    po = os.path.splitext(p)[0] + ".pdf"
    print(po)
    Plot(u, po)
