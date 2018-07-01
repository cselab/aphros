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

from matplotlib.colors import LinearSegmentedColormap

# u: 2d numpy array
# po: output path
def PlotFieldGray(ax, u, vmin=None, vmax=None):
    cmap = LinearSegmentedColormap.from_list('mycmap', ['#c8c5bd', '#787672'])
    ax.imshow(np.flipud(u), extent=(0, 1, 0, 1), interpolation='nearest',
              vmin=vmin, vmax=vmax, cmap=cmap)

# u: 2d numpy array
# po: output path
def PlotField(ax, u, vmin=None, vmax=None):
    ax.imshow(np.flipud(u), extent=(0, 1, 0, 1), interpolation='nearest',
              vmin=vmin, vmax=vmax)

# u: 2d numpy array
# po: output path
def PlotFieldBwr(ax, u, vmin=None, vmax=None):
    ax.imshow(np.flipud(u), extent=(0, 1, 0, 1), interpolation='nearest',
              vmin=vmin, vmax=vmax, cmap=plt.get_cmap("bwr"))

# sx, sy: number of cells
def PlotGrid(ax, x1, y1):
    ax.set_xticks(x1)
    ax.set_yticks(y1)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    ax.grid(True, lw=0.5, c='0.', alpha=0.08)

def PlotLines(ax, xa, ya, xb, yb):
    xa = xa.flatten()
    ya = ya.flatten()
    xb = xb.flatten()
    yb = yb.flatten()
    xy = np.vstack((xa, xb, ya, yb))
    xy = xy.T
    xy = xy.reshape((xa.size * 2, 2))
    ax.plot(*xy, c='k', alpha=0.9, lw=1)

# assume 0 < nx < ny, a < 0
def GetLineEnds(nx, ny, a):
    xl = [(a + 0.5 * ny) / nx, (a + 0.5 * nx) / ny]
    xr = [(a - 0.5 * ny) / nx, (a - 0.5 * nx) / ny]

    e = []
    if -0.5 <= xl[0] and xl[0] <= 0.5:
        e.append([xl[0], -0.5])
    if -0.5 <= xr[0] and xr[0] <= 0.5:
        e.append([xr[0], 0.5])
    if len(e) < 2 and -0.5 <= xl[1] and xl[1] <= 0.5:
        e.append([-0.5, xl[1]])
    if len(e) < 2 and -0.5 <= xr[1] and xr[1] <= 0.5:
        e.append([0.5, xr[1]])

    if len(e) == 1: # if only one point found, set second to the same
        e.append(e[0])

    return e


# xc,yc: cell centers
# a: line constants
# nx, ny: unit normals
# hx, hy: cell size
# Returns:
# xa, ya, xb, yb: line end points (a,b) on cell edges
# Equation of line:
# (x-xc)/h dot n = a
def GetLines(xc, yc, a, nx, ny, hx, hy, u):
    xc = xc.flatten()
    yc = yc.flatten()
    nx = nx.flatten()
    ny = ny.flatten()
    a = a.flatten()
    u = u.flatten()

    xa = []
    ya = []
    xb = []
    yb = []
    for i in range(len(xc)):
        th = 1e-5
        if u[i] > th and u[i] < 1. - th:
            e = GetLineEnds(nx[i] * hx, ny[i] * hy, a[i])
            if len(e) == 2:
                xa.append(xc[i] + e[0][0] * hx)
                ya.append(yc[i] + e[0][1] * hy)
                xb.append(xc[i] + e[1][0] * hx)
                yb.append(yc[i] + e[1][1] * hy)

    xa = np.array(xa)
    ya = np.array(ya)
    xb = np.array(xb)
    yb = np.array(yb)

    return xa, ya, xb, yb

# kk: list of fields
# ll: list of labels
def CmpCurv(x, y, u, kk, ll, po):
    global reff, hx

    cx = (x * u).sum() / u.sum()
    cy = (y * u).sum() / u.sum()
    th = 1e-3
    ii = np.where((u > th) & (u < 1. - th))
    x = x[ii];   y = y[ii]
    dx = x - cx; dy = y - cy
    an = np.arctan2(dy, dx)
    deg = np.degrees(an)

    s = np.argsort(an)


    # exact curvature
    ane = np.linspace(-np.pi, np.pi, 200)
    dege = np.degrees(ane)
    cx,cy,cz,rx,ry = np.loadtxt('b.dat')
    ke = rx*ry / ((ry * np.cos(ane)) ** 2 + (rx * np.sin(ane)) ** 2) ** (3. / 2.)

    #kea = 1. / reff  # average curvature from area
    kea = ke.mean()

    fig, ax = plt.subplots()
    #ax.axhline(y=1., c="0.5")
    ax.plot(dege, ke / kea, label="exact", c="0.5")

    for k,l in zip(kk,ll):
        if k is not None:
            print("plot l={:}, median={:}".format(l, np.max(k)))
            ax.plot(deg[s], k[ii][s] / kea, label=l)

    ax.set_xlabel(r"angle [deg]")
    ax.set_ylabel(r"normalized curvature")
    q=180.
    ax.set_xticks(np.arange(-q, q*1.1, 90.))
    ax.set_xlim(-q, q)
    ax.legend()
    ax.set_ylim(0., 2.)
    ax.grid()
    #plt.title("R/h={:.3f}".format(reff / hx))
    plt.title("rx/h={:.3f}, ry/h={:.3f}".format(rx / hx, ry / hx))
    fig.tight_layout()
    fig.savefig(po, dpi=300)
    plt.close()

# Plot particle strings
# ax: axes to plot on
# p: path to csv with columns x,y,z,c (c: cell index)
# sk: plot every sk string
# Output:
# plots particles on ax with color c
def PlotPart(ax, p, sk=1):
    d = np.loadtxt(p, skiprows=1, delimiter=',')
    if d.size == 0:
        return None
    x,y,z,c = d.T
    c = c.astype(int)


    # separate particles strings by cell
    tt = dict()  # lists of indices in x
    for i in range(len(c)):
        ci = c[i]
        if not ci in tt:
            tt[ci] = []
        tt[ci].append(i)

    # map cell index to color
    nc = 16
    c = (c * (2 ** 31 - 1) % nc).astype(float) / (nc - 1)

    cmap = plt.get_cmap("Set1")

    # plot strings
    n = 0
    for i in tt:
        if i % sk == 0:
            ti = np.array(tt[i])
            cl = cmap(c[ti[0]])
            # connecting lines
            ax.plot(x[ti], y[ti], c=cl, zorder=10, lw=1, alpha=0.5)
            # points
            ax.scatter(x[ti], y[ti], c=cl, s=2, lw=0.3, zorder=11, edgecolor='black')
        n += 1


pre = 'u'
pp = sorted(glob.glob(pre + "*.dat"))
#pp = pp[1:2]
for p in pp:
    suf = re.findall(pre + "(.*)", p)[0]
    print("p={:}, suf={:}".format(p, suf))

    u = Get2d(Read(p))

    [sx, sy] = u.shape # cells

    # check if gerris (odd sx)
    ge = (sx % 2 == 1)

    if ge:
        sx -= 1
        sy -= 1

    hx = 1. / sx
    hy = 1. / sy

    # equivalent radius
    reff = ((u.sum() * hx * hy) / np.pi) ** 0.5
    print("reff={:}".format(reff))

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


    # invert curvature if gerris
    if ge:
        k = -k

    # volume fraction u
    po = os.path.splitext(p)[0] + ".pdf"
    print(po)
    fig, ax = PlotInit()
    PlotGrid(ax, xn1, yn1)
    vmax = abs(u).max()
    PlotFieldGray(ax, u, vmin=0., vmax=1.)
    #PlotField(ax, np.clip(u, 0., 1.))
    if all([e is not None for e in [a, nx, ny]]):
        l = GetLines(x, y, a, nx, ny, hx, hy, u)
        PlotLines(ax, *l)
    # plot partilces
    n = re.findall("_(\d*)\.", suf)
    if n:
        pa = "partit.{:}.csv".format(n[0])
        if os.path.isfile(pa):
            print(pa)
            PlotPart(ax, pa, sk=5)
    # save
    PlotSave(fig, ax, po)

    # curvature k
    if k is not None:
        po = os.path.splitext('k' + suf)[0] + ".pdf"
        print(po)
        fig, ax = PlotInit()
        PlotGrid(ax, xn1, yn1)

        # average curvature
        ke = 1. / reff  # exact

        # exclude points with zero curvature (for ch)
        k[np.where(k == 0.)] = np.nan

        PlotFieldBwr(ax, k / ke, vmin=0.8, vmax=1.2)
        # plot lines
        if all([e is not None for e in [a, nx, ny]]):
            l = GetLines(x, y, a, nx, ny, hx, hy, u)
            PlotLines(ax, *l)

        PlotSave(fig, ax, po)

    # plot curvature comparison
    kh = Get2d(Read('kh' + suf))
    kp = Get2d(Read('kp' + suf))
    kg = Get2d(Read('gerris_k/k0011.dat'))
    if kg is not None:
        kg = -kg
    if kh is not None and kp is not None:
        po = os.path.splitext('kcmp' + suf)[0] + ".pdf"
        kk = [kh, kp, kg]
        ll = ["height", "particles", "gerris"]
        CmpCurv(x, y, u, kk, ll, po)

