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
    plt.close()
    fig, ax = plt.subplots(figsize=(5,5))
    return fig, ax

def PlotSave(fig, ax, po):
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(po, dpi=300)
    plt.close()

# u: 2d numpy array
# po: output path
def PlotField(ax, u, vmin=None, vmax=None):
    ax.imshow(np.flipud(u), extent=(0, 1, 0, 1), interpolation='nearest',
              vmin=vmin, vmax=vmax)

# sx, sy: number of cells
def PlotGrid(ax, x1, y1):
    ax.set_xticks(x1)
    ax.set_yticks(y1)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    ax.grid(True)

def PlotLines(ax, xa, ya, xb, yb):
    xa = xa.flatten()
    ya = ya.flatten()
    xb = xb.flatten()
    yb = yb.flatten()
    xy = np.vstack((xa, xb, ya, yb))
    xy = xy.T
    xy = xy.reshape((xa.size * 2, 2))
    ax.plot(*xy, c='r')

# xc,yc: cell centers
# a: line constants
# nx, ny: unit normals
# hx, hy: cell size
# Returns:
# xa, ya, xb, yb: line end points (a,b) on cell edges
# Equation of line:
# (x-xc)/h dot n = a
def GetLines(xc, yc, a, nx, ny, hx, hy, u):
    def ClipX(x, y, tx, ty, xmin, xmax):
        xp = np.clip(x, xmin, xmax)
        dx = xp - x
        tx = np.where(tx == 0., 1e-10, tx)
        dy = dx / tx * ty
        return x + dx, y + dy
    def ClipY(x, y, tx, ty, ymin, ymax):
        y, x = ClipX(y, x, ty, tx, ymin, ymax)
        return x, y
    def Clip(x, y, tx, ty, xmin, xmax, ymin, ymax):
        x, y = ClipX(x, y, tx, ty, xmin, xmax)
        x, y = ClipY(x, y, tx, ty, ymin, ymax)
        return x, y
    xc = xc.flatten()
    yc = yc.flatten()
    nx = nx.flatten()
    ny = ny.flatten()
    a = a.flatten()
    u = u.flatten()
    # tangent
    tx = ny
    ty = -nx
    # line center
    xlc = xc + a * nx
    ylc = yc + a * ny
    # end points
    xa = xlc - tx * hx
    ya = ylc - ty * hy
    xb = xlc + tx * hx
    yb = ylc + ty * hy
    # cell bounds
    xcm = xc - hx * 0.5
    ycm = yc - hy * 0.5
    xcp = xc + hx * 0.5
    ycp = yc + hy * 0.5
    # clip
    xa, ya = Clip(xa, ya, tx, ty, xcm, xcp, ycm, ycp)
    xb, yb = Clip(xb, yb, tx, ty, xcm, xcp, ycm, ycp)

    th = 1e-5
    w = np.where((u > th) & (u < 1. - th))[0]
    xa = xa[w]
    ya = ya[w]
    xb = xb[w]
    yb = yb[w]

    return xa, ya, xb, yb

def CmpCurv(x, y, u, k, kp, po):
    cx = (x * u).sum() / u.sum()
    cy = (y * u).sum() / u.sum()
    th = 1e-3
    ii = np.where((u > th) & (u < 1. - th))
    x = x[ii];   y = y[ii]
    dx = x - cx; dy = y - cy
    an = np.arctan2(dy, dx)
    deg = np.degrees(an)
    k = k[ii]
    kp = kp[ii]

    s = np.argsort(an)

    [sx, sy] = u.shape
    hx = 1. / sx
    hy = 1. / sy
    r = ((u.sum() * hx * hy) / np.pi) ** 0.5
    ke = 1. / r

    fig, ax = plt.subplots()
    ax.plot(deg[s], k[s] / ke, label="k")
    ax.plot(deg[s], kp[s] / ke, label="kp")
    ax.legend()
    fig.tight_layout()
    fig.savefig(po, dpi=300)
    plt.close()



pre = 'u'
pp = sorted(glob.glob(pre + "*.dat"))
#pp = pp[1:2]
for p in pp[10:11]:
    print(p)
    suf = re.findall(pre + "(.*)", p)[0]
    u = Get2d(Read(p))
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
    a = Get2d(Read('a' + suf)) # alpha
    nx = Get2d(Read('nx' + suf)) # normal
    ny = Get2d(Read('ny' + suf))
    k = Get2d(Read('k' + suf))

    # volume fraction u
    po = os.path.splitext(p)[0] + ".pdf"
    print(po)
    fig, ax = PlotInit()
    PlotGrid(ax, xn1, yn1)
    vmax = abs(u).max()
    PlotField(ax, u, vmin=-vmax, vmax=vmax)
    #PlotField(ax, np.clip(u, 0., 1.))
    if all([e is not None for e in [a, nx, ny]]):
        l = GetLines(x, y, a, nx, ny, hx, hy, u)
        PlotLines(ax, *l)
    PlotSave(fig, ax, po)

    # curvature k
    if k is not None:
        po = os.path.splitext('k' + suf)[0] + ".pdf"
        print(po)
        fig, ax = PlotInit()
        PlotGrid(ax, xn1, yn1)
        vmax = 1. / 0.2 * 2.
        PlotField(ax, k, vmin=-vmax, vmax=vmax)
        if all([e is not None for e in [a, nx, ny]]):
            l = GetLines(x, y, a, nx, ny, hx, hy, u)
            PlotLines(ax, *l)
        PlotSave(fig, ax, po)

    kp = Get2d(Read('kp' + suf))
    if k is not None and kp is not None:
        po = os.path.splitext('kcmp' + suf)[0] + ".pdf"
        CmpCurv(x, y, u, k, kp, po)
