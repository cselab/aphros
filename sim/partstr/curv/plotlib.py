#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re
from matplotlib.colors import LinearSegmentedColormap

# natural sort
def natkey(s, _nsre=re.compile('([0-9]+)')):
      return [int(text) if text.isdigit() else text.lower()
                      for text in re.split(_nsre, s)]

def natsorted(v):
  return sorted(v, key=natkey)

def Glob(d, fld):
    return natsorted(glob.glob(os.path.join(d, "{:}*.dat".format(fld))))

# Read uniform grid data
# p: path
# Format:
# <Nx> <Ny> <Nz>
# <u[0,0,0]> <u[0,0,1]> ...
# Return:
# array of shape (Nx, Ny, Nz) (reversed dimensions)
# None if file not found
def ReadArray(p):
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
        # data x,y,z
        u = np.transpose(u)
        return u

# Converts field to image
def GetImg(u):
    s = u.shape
    assert len(s) in [2, 3]
    if len(s) == 3:
        u = u[:,:,0]
    return np.flipud(u.T)

def PlotInitSq():
    fig, ax = plt.subplots(figsize=(3,3))
    ax.set_aspect('equal')
    return fig, ax

def PlotInit():
    fig, ax = plt.subplots(figsize=(4,3))
    return fig, ax

def PlotInit3():
    fig = plt.figure(figsize=(4,3))
    ax = fig.gca(projection='3d')
    return fig, ax


def PlotSave(fig, ax, po):
    fig.tight_layout()
    fig.savefig(po, dpi=300)
    plt.close()

from matplotlib.colors import LinearSegmentedColormap

# u: 2d numpy array
# po: output path
def PlotFieldGray(ax, u, vmin=None, vmax=None):
    cmap = LinearSegmentedColormap.from_list('mycmap', ['#c8c5bd', '#787672'])
    ax.imshow(GetImg(u), extent=(0, 1, 0, 1), interpolation='nearest',
              vmin=vmin, vmax=vmax, cmap=cmap)

# u: 2d numpy array
# po: output path
def PlotField(ax, u, vmin=None, vmax=None):
    ax.imshow(GetImg(u), extent=(0, 1, 0, 1), interpolation='nearest',
              vmin=vmin, vmax=vmax)

# u: 2d numpy array
# po: output path
def PlotFieldBwr(ax, u, vmin=None, vmax=None):
    ax.imshow(GetImg(u), extent=(0, 1, 0, 1), interpolation='nearest',
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
        th = 1e-10
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
    x,y,z,c = d.T[0:4,:]
    c = c.astype(int)


    # separate particles strings by cell
    tt = dict()  # lists of indices in x, cell is key
    for i in range(len(c)):
        ci = c[i]
        if not ci in tt:
            tt[ci] = []
        tt[ci].append(i)

    # map cell index to [0:nc]
    nc = 16
    #c = (c * (2 ** 31 - 1) % nc).astype(float) / (nc - 1)
    cu = np.unique(c)
    np.random.seed(0)
    rnd = np.random.choice(np.arange(nc), len(cu))
    mcl = dict(zip(cu, rnd))

    cmap = plt.get_cmap("Set1")

    # plot strings
    for i in tt:
        if np.random.randint(sk) == 0:
            ti = np.array(tt[i])
            cl = cmap(mcl[i])
            # connecting lines
            ax.plot(x[ti], y[ti], c=cl, zorder=10, lw=1, alpha=0.5)
            # points
            ax.scatter(x[ti], y[ti], c=cl, s=2, lw=0.3, zorder=11, edgecolor='black')
            #ax.scatter(x[ti], y[ti], c=cl, s=0.5, lw=0.2, zorder=11, alpha=0.8, edgecolor='black') # XXX

def IsGerris(s):
    # check if gerris (odd nx)
    return s[0] % 2 == 1

# s: shape of data array
# Returns:
# x1,y1,z1: coordinates of cell centers
# hx,hy,hz: mesh steps (hz=1 if 2d))
def GetGeom(s):
    [nx, ny, nz] = s
    # cells
    ge = IsGerris(s)
    sx,sy,sz = (nx,ny,nz) if not ge else (nx-1,ny-1,nz-1)
    hx = 1. / sx
    hy = 1. / sy
    hz = 1. / sz if sz > 0 else 1.
    # mesh
    if ge:
        x1 = (0. + np.arange(nx)) * hx
        y1 = (0. + np.arange(ny)) * hy
        z1 = (0. + np.arange(nz)) * hz
    else:
        x1 = (0.5 + np.arange(nx)) * hx
        y1 = (0.5 + np.arange(ny)) * hy
        z1 = (0.5 + np.arange(nz)) * hz
    return x1,y1,z1,hx,hy,hz

def GetMesh(x, y, z):
    return np.meshgrid(x, y, z, indexing='ij')

# assuming uniform mesh
def GetMeshStep(s):
    if IsGerris(s):
        return 1. / (s[0] - 1)
    return 1. / s[0]

def GetDim(s):
    if len(s) == 2 or s[2] == 1:
        return 2
    return 3

# Returns mesh nodes
# assume [0,1]x[0,1]x[0,1] domain
def GetMeshNodes(hx, hy, hz):
    xn1 = np.arange(0., 1. + hx * 0.5, hx)
    yn1 = np.arange(0., 1. + hy * 0.5, hy)
    zn1 = np.arange(0., 1. + hz * 0.5, hz)
    return xn1, yn1, zn1

# Radius of circle form area
def GetReff2(a):
    return (a / np.pi) ** 0.5

# Radius of sphere from volume
def GetReff3(v):
    return (v / np.pi * 3. / 4.) ** (1. / 3.)

# Returns effective radius
def GetReff(vf, dim):
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    # cell volume
    vc = hx * hy * hz
    # integral of vf
    v  = vf.sum() * vc
    return GetReff2(v / hz) if dim == 2 else GetReff3(v)

# Returns path template from sample path
def GetPathTemplate(p):
    d = os.path.dirname(p)
    b = os.path.basename(p)
    b = re.sub("[^\d]*", "{:}", b, count=1)
    b = os.path.splitext(b)[0] + ".{:}"
    pt = os.path.join(d, b)
    return pt

# Returns full path to field
# pt: path template, example: "ch/{:}_0.dat"
# fld: field name
def GetFieldPath(pt, fld, ext="dat"):
    if fld == "partit":
        return pt.format(fld + "_", "csv")
    return pt.format(fld + "_", ext)

def ReadField(pt, fld):
    p = GetFieldPath(pt, fld)
    return ReadArray(GetFieldPath(pt, fld))

# Save figure with title only
# t: title
# po: output path
def FigTitle(t, po):
    fig, ax = plt.subplots(figsize=(5,0.5))
    ax.text(0., 0., t, fontsize=15)
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    #fig.tight_layout()
    #fig.savefig(po, bbox_inches='tight')
    fig.savefig(po)
    plt.close()

# Computes norms: max, L1, L2
# u: array
# Returns:
# m, l1, l2
def GetNorm(k):
    k = abs(k)
    m = k.max()
    l1 = k.mean()
    l2 = (k ** 2).mean() ** 0.5
    return m, l1, l2

def LoadBub():
    return np.loadtxt("b.dat")

# Returns mean curvature of surface z=h(x,y) at point (x,y,z)
# h: function h(x,y)
# x,y: coordinates
# dx: step to compute derivative
def GetCurv(h, x, y, dx=1e-3):
    dy = dx
    # returns df/dx of function f(x,y)
    def Dx(f):
        return lambda x,y: (f(x + dx, y) - f(x - dx, y)) / (2. * dx)
    # returns df/dy of function f(x,y)
    def Dy(f):
        return lambda x,y: (f(x, y + dy) - f(x, y - dy)) / (2. * dy)

    hx = Dx(h)(x,y)
    hy = Dy(h)(x,y)
    hxx = Dx(Dx(h))(x,y)
    hxy = 0.5 * (Dx(Dy(h))(x,y) + Dy(Dx(h))(x,y))
    hyy = Dy(Dy(h))(x,y)

    #a = (1. + hx ** 2) * hyy - 2. * hx * hy * hxy + (1. + hy ** 2) * hxx
    #c = (hx ** 2 + hy ** 2 + 1.) ** (3. / 2)
    #return -a / c

    a = hx ** 2 * hxx + 2. * hx * hy * hxy + hy ** 2 * hyy
    b = (hx ** 2 + hy ** 2 + 1.) * (hxx + hyy)
    c = (hx ** 2 + hy ** 2 + 1.) ** (3. / 2)
    return (a - b) / c

def P(x, lbl):
    print("{:}: shape={:} min={:}, max={:}, avg={:}".format(
        lbl, x.shape, x.min(), x.max(), x.mean()))

# Returns curvature of ellipse
# (x/rx)**2 + (y/ry)**2 = 1 at point x
def GetEllip2Curv0(x, rx, ry):
    dx = rx * 1e-3
    x = np.clip(x, -rx + dx * 2, rx - dx * 2) # XXX clip
    return GetCurv(
        lambda xx,zz:
        max(0., 1. - (xx / rx) ** 2) ** 0.5 * ry,
        x, 0., dx=dx)

# Returns curvature of ellipse
# (x/rx)**2 + (y/ry)**2 = 1 at point x,y
def GetEllip2Curv(x, y, rx, ry):
    if abs(x) < abs(y):
        return GetEllip2Curv0(x, rx, ry)
    return GetEllip2Curv0(y, ry, rx)

# Returns curvature of ellipsoid
# (x/rx)**2 + (y/ry)**2 + (z/rz) ** 2 = 1 at point x,y
def GetEllip3Curv0(x, y, rx, ry, rz):
    dd = 1e-3
    dx = rx * dd
    dy = ry * dd
    return GetCurv(
        lambda xx,yy: max(0., 1. - (xx / rx) ** 2 - (yy / ry) ** 2) ** 0.5 * rz,
        x, y, dx=dx)

# Returns curvature of ellipsoid
# (x/rx)**2 + (y/ry)**2 + (z/rz) ** 2 = 1 at point x,y,z
def GetEllip3Curv(x, y, z, rx, ry, rz):
    if abs(z) >= abs(x) and abs(z) >= abs(y): # z(x,y)
        return GetEllip3Curv0(x, y, rx, ry, rz)
    elif abs(y) >= abs(x) and abs(y) >= abs(z): # y(x,z)
        return GetEllip3Curv0(x, z, rx, rz, ry)
    # x(y,z)
    return GetEllip3Curv0(y, z, ry, rz, rx)

# Evaluates exact curvature of ellipsoid
# dim: dimension, 2 or 3
# x,y,z: points
# Returns:
# k: curvature of shape x
def GetExactK(dim, x, y, z):
    cx,cy,cz,rx,ry,rz = LoadBub()

    dx = x - cx; dy = y - cy; dz = z - cz

    k = np.zeros_like(x)

    if dim == 2:
        u = np.arctan2(dy, dx)
        r = ((np.cos(u) / rx) ** 2 + (np.sin(u) / ry) ** 2) ** (-0.5)
        dx = r * np.cos(u)
        dy = r * np.sin(u)
        for i,wx,wy in zip(range(len(dx)),dx,dy):
            k[i] = GetEllip2Curv(wx, wy, rx, ry)

    if dim == 3:
        u = np.arctan2(dy, dx)
        v = np.arctan2(dz, (dx ** 2 + dy ** 2) ** 0.5)
        r = ((np.cos(u) * np.cos(v) / rx) ** 2 +
             (np.sin(u) * np.cos(v) / ry) ** 2 +
             (np.sin(v) / rz) ** 2) ** (-0.5)
        dx = r * np.cos(u) * np.cos(v)
        dy = r * np.sin(u) * np.cos(v)
        dz = r * np.sin(v)
        for i,wx,wy,wz in zip(range(len(dx)),dx,dy,dz):
            k[i] = GetEllip3Curv(wx, wy, wz, rx, ry, rz)

    return k


# Histogram of curvature
def FigHistK(vf, kk, ll, po, title=None):
    dim = GetDim(vf.shape)
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    x,y,z = GetMesh(x1, y1, z1)

    # interface cells
    th = 1e-8
    ii = np.where((vf > th) & (vf < 1. - th))
    x = x[ii]; y = y[ii]; z = z[ii]

    # exact curvature
    ke = GetExactK(dim, x, y, z)
    # average curvature
    kea = ke.mean()

    fig, ax = PlotInit()
    for k,l in zip(kk,ll):
        if IsGerris(vf.shape):
            k *= -1
        k = k[ii]
        er = (k - ke) / kea  # error
        h,b = np.histogram(er, 200, range=(-1., 1.), density=False)
        bc = (b[1:] + b[:-1]) * 0.5
        ax.plot(bc, h, label=l)

    ax.axvline(x=0., label="exact", c="0.5", ls='--')
    ax.set_ylim(0)
    ax.set_xlabel(r"curvature error")
    ax.legend()
    if title is not None:
        ax.set_title(title)
    ax.grid()
    PlotSave(fig, ax, po)

    # write error: label rx ry max l1 l2
    cx,cy,cz,rx,ry,rz = LoadBub()
    with open("er", 'w') as f:
        for k,l in zip(kk,ll):
            if k is not None:
                k = k[ii]
                m,l1,l2 = GetNorm((k - ke) / kea)
                f.write("{:} {:} {:} {:} {:} {:}\n".format(
                        l, rx/hx, ry/hx, m, l1, l2))
