#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re

# Read uniform grid data
# p: path
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

# u: numpy array (2d or 3d with shape[2]==1)
def Get2d(u):
    if len(u.shape) == 2:
        return u
    else:
        s = u.shape
        assert len(s) == 3 and s[2] == 1
        return u.reshape((s[0], s[1]))

def PlotInit():
    fig, ax = plt.subplots(figsize=(5,5))
    return fig, ax

def PlotSave(fig, ax, po):
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(po, dpi=300)
    plt.close()

# u: 2d numpy array
# po: output path
def PlotField(ax, u):
    u = np.clip(u, 0., 1.)
    ax.imshow(np.flipud(u), extent=(0, 1, 0, 1), interpolation='nearest')

# sx, sy: number of cells
def PlotGrid(ax, x1, y1):
    ax.set_xticks(x1)
    ax.set_yticks(y1)
    #ax.xaxis.set_ticklabels([])
    #ax.yaxis.set_ticklabels([])
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    ax.grid(True)

def PlotLines(ax, xa, xb, ya, yb):
    ax.plot(xa, xb, ya, yb)


pre = 'u'
pp = sorted(glob.glob(pre + "*.dat"))
pp = pp[1:2]
for p in pp:
    suf = re.findall(pre + "(.*)", p)[0]
    u = Read(p)
    u = Get2d(u)
    [sx, sy] = u.shape
    hx = 1. / sx
    hy = 1. / sy
    # nodes, 1 means 1d
    xn1 = np.arange(sx + 1) * hx
    yn1 = np.arange(sy + 1) * hy
    # cell centers
    x1 = (0.5 + np.arange(sx)) * hx
    y1 = (0.5 + np.arange(sy)) * hy
    x, y = np.meshgrid(x1, y1)
    a = Read('a' + suf) # alpha

    po = os.path.splitext(p)[0] + ".pdf"
    print(po)
    fig, ax = PlotInit()
    PlotGrid(ax, xn1, yn1)
    PlotField(ax, u)
    #PlotLines(ax, x, y, x, y)
    PlotSave(fig, ax, po)
