#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re
import stream

# Read uniform grid data
# p: path
# Format:
# <Nx> <Ny> <Nz>
# <u[0,0,0]> <u[0,0,1]> ...
# Return:
# array of shape (Nz, Ny, Nx)
# None if file not found
def Read(p):
    if not os.path.isfile(p):
        return None
    with open(p) as f:
        ll = f.readlines()
        # shape x,y,z
        s = np.array(ll[0].split(), dtype=int)
        # shape z,y,x
        ss = tuple(reversed(s))
        # data flat
        u = np.array(ll[1].split(), dtype=float)
        # data z,y,x
        u = u.reshape(ss)
        return u

# u: numpy array (2d or 3d slice)
def Get2d(u):
    if u is None:
        return None
    s = u.shape
    if len(s) == 2:
        return u
    else:
        assert len(s) == 3
        return u[0,:,:].reshape((s[1], s[2]))

def PlotInit():
    fig, ax = plt.subplots(figsize=(5,5))
    return fig, ax

def PlotSave2(fig, ax, po):
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(po, dpi=300)
    plt.close()

def PlotSave1(fig, ax, po):
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
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    ax.grid(True)

def PlotStream(x1, y1, vx, vy, o):
    psi = stream.stream(vx, vy)
    fig, ax = PlotInit()
    plt.contour(x1, y1, psi, 36, colors='k',
                linestyles="solid", linewidths=1.)
    ax.imshow(np.flipud(p), extent=(0, 1, 0, 1), interpolation='nearest')
    PlotSave2(fig, ax, o)

def PlotVy(ax, d, l=None):
    if (l is None):
        l = d

    os.path.isdir(d) or os.mkdir(d)
    d = os.path.join(d, "")
    x = np.loadtxt(d + "x")
    vy = np.loadtxt(d + "vy")
    ax.plot(x, vy, label=l)

def SaveVy(x1, vy, d):
    os.path.isdir(d) or os.mkdir(d)
    d = os.path.join(d, "")
    ny = vx.shape[1];
    vy1 = vy[ny // 2,:]
    np.savetxt(d + "x", x1)
    np.savetxt(d + "vy", vy1)

def FigVy(x, vy, o):
    SaveVy(x, vy, "ch")

    plt.close()
    fig, ax = PlotInit()

    PlotVy(ax, "ref")
    PlotVy(ax, "ch64q")
    PlotVy(ax, "ch", "ch,nx={:}".format(x.size))

    ax.grid()
    ax.legend(loc="best")
    PlotSave1(fig, ax, o)

pre = 'p'
ff = sorted(glob.glob(pre + "*.dat"))[0:]
#pp = pp[1:2]
for f in ff:
    suf = re.findall(pre + "(.*)", f)[0]
    vx = Get2d(Read('vx' + suf))
    vy = Get2d(Read('vy' + suf))
    p = Get2d(Read('p' + suf))
    [sx, sy] = vx.shape
    hx = 1. / sx
    hy = 1. / sy
    # nodes, 1 means 1d
    xn1 = np.arange(sx + 1) * hx
    yn1 = np.arange(sy + 1) * hy
    # cell centers
    x1 = (0.5 + np.arange(sx)) * hx
    y1 = (0.5 + np.arange(sy)) * hy
    x, y = np.meshgrid(x1, y1)

    po = os.path.splitext(f)[0] + ".pdf"
    print(po)
    PlotStream(x1, y1, vx, vy, po)

    po = os.path.splitext(f)[0] + "_vy.pdf"
    print(po)
    FigVy(x1, vy, po)
