#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re
import scipy.interpolate
import sys

# natural sort
def natkey(s, _nsre=re.compile('([0-9]+)')):
      return [int(text) if text.isdigit() else text.lower()
                      for text in re.split(_nsre, s)]

def natsorted(v):
  return sorted(v, key=natkey)

def GetLast(p):
    return natsorted(glob.glob(p))[-1]


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

def Read2d(p):
    return Get2d(Read(p))

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

# vvx: list of velocity magnitude
# vvy: list of velocity magnitude
# ll: list of labels
def Cmp(x, y, x1, y1, vf, vvx, vvy, ll, po, vx0, vy0):
    global reff, hx, ge

    cx = (x * vf).sum() / vf.sum()
    cy = (y * vf).sum() / vf.sum()

    rmax = min(cx, 1. - cx, cy, 1. - cy) * 0.9
    r1 = np.linspace(0., rmax, 100)  # radius
    t1 = np.linspace(0., 2. * np.pi, 200)  # angle
    # mesh of radius and angle, index: [r,t]
    tt,rr = np.meshgrid(t1, r1)

    xx = cx + rr * np.cos(tt)
    yy = cy + rr * np.sin(tt)

    fig, ax = plt.subplots()
    #ax.axhline(y=0., c='k')
    ax.axvline(x=reff, c='k', ls='--')

    for vx,vy,l in zip(vvx,vvy,ll):
        if vx is not None and vy is not None:
            # velocity magnitude
            m = np.sqrt((vx - vx0) ** 2 + (vy - vy0) ** 2)
            m = m[:len(x1),:len(y1)]
            print("plot l={:}, median={:}".format(l, np.median(m)))
            mrt = scipy.interpolate.interpn((x1, y1), m, (xx, yy),
                                          bounds_error=False, fill_value=0.)
            mr = np.mean(mrt, axis=1)
            ax.plot(r1, mr, label=l)

    ax.set_xlabel(r"distance from bubble center")
    ax.set_ylabel(r"radial mean of velocity magnitude")
    #ax.set_xticks(np.arange(-q, q*1.1, 90.))
    ax.legend()
    ax.set_ylim(0.)
    ax.set_xlim(0., rmax)
    ax.grid()
    plt.title("R={:.3f}, R/h={:.3f}".format(reff, reff / hx))
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

# s: shape of data array
# Returns:
# x1,y1: coordinates of cell centers
def GetGeom(s):
    [nx, ny] = s
    # check if gerris (odd nx)
    ge = (nx % 2 == 1)
    # cells
    sx,sy = (nx,ny) if not ge else (nx-1,ny-1)
    hx = 1. / sx
    hy = 1. / sy
    # mesh
    if ge:
        x1 = (0. + np.arange(nx)) * hx
        y1 = (0. + np.arange(nx)) * hy
    else:
        x1 = (0.5 + np.arange(nx)) * hx
        y1 = (0.5 + np.arange(nx)) * hy
    #x, y = np.meshgrid(x1, y1)
    return x1,y1,hx,hy

def GetMeshNodes(s):
    x1,y1,hx,hy = GetGeom(s)
    xn1 = np.arange(0., 1. + hx * 0.5, hx)
    yn1 = np.arange(0., 1. + hx * 0.5, hy)
    return xn1, yn1

# pre: prefix
# d: directory
def GetFiles(d, pre):
    pp = natsorted(glob.glob("{:}/{:}*.dat".format(d, pre)))
    return pp

# Extracts trajectory of center of mass.
# pp: list of paths to volume fraction fields
# Returns:
# x,y: arrays
def GetTraj(pp):
    # result
    cx = [] ; cy = []
    for i,p in enumerate(pp):
        # report
        sys.stdout.write("{:} ".format(i))
        sys.stdout.flush()
        # read array
        vf = Read2d(p)
        x1,y1,hx,hy = GetGeom(vf.shape)
        x,y = np.meshgrid(x1, y1)
        # compute center
        cx0 = (x * vf).sum() / vf.sum()
        cy0 = (y * vf).sum() / vf.sum()
        cx.append(cx0)
        cy.append(cy0)

    cx = np.array(cx)
    cy = np.array(cy)
    print()
    return cx,cy

# Plots trajectories
# xx,yy: list of arrays
# ll: labels
# s: shape of of data array
def PlotTraj(xx, yy, ll, s):
    global reff, hx
    fig, ax = PlotInit()
    ax.set_aspect('equal')
    for x,y,l in zip(xx, yy, ll):
        ax.plot(x, y, label=l)
    #ax.set_xlabel(r"x")
    #ax.set_ylabel(r"y")
    ax.legend()
    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    xn1,yn1 = GetMeshNodes(s)
    PlotGrid(ax, xn1, yn1)
    plt.title("R/h={:.3f}".format(reff / hx))
    po = 'cmp.pdf'
    PlotSave(fig, ax, po)

# Returns effective radius
def GetReff(vf):
    x1,y1,hx,hy = GetGeom(vf.shape)
    return ((vf.sum() * hx * hy) / np.pi) ** 0.5

# Exact trajectory.
# x0,y0: initial center
# vx0,vy0: velocity
# t: time
# n: frames
def GetTrajE(x0, y0, vx0, vy0, t, n=100):
    x = np.linspace(x0, x0 + vx0 * t, n)
    y = np.linspace(y0, y0 + vy0 * t, n)
    return x,y

dd = ["./", "gerris"]
ll = ["ch", "gerris"]
xx = []
yy = []

global reff, hx

pre = "vf"
for d,l in zip(dd, ll):
    pp = GetFiles(d, pre)
    vf = Read2d(pp[0])
    x1,y1,hx,hy = GetGeom(vf.shape)
    reff = GetReff(vf)
    x,y = GetTraj(pp)
    xx.append(x)
    yy.append(y)

vf = Read2d(GetFiles(dd[0], pre)[0])
x1,y1,hx,hy = GetGeom(vf.shape)
reff = GetReff(vf)

# exact trajectory
x,y = GetTrajE(0.3, 0.3, 0.4, 0.3, 1.)
xx = [x] + xx
yy = [y] + yy
ll = ["exact"] + ll

# Plot trajectories
PlotTraj(xx, yy, ll, vf.shape)

# Plot volume fraction
xn1,yn1 = GetMeshNodes(vf.shape)
po = "vf.pdf"
print(po)
fig, ax = PlotInit()
PlotGrid(ax, xn1, yn1)
PlotFieldGray(ax, vf, vmin=0., vmax=1.)
PlotSave(fig, ax, po)

