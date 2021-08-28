#!/usr/bin/env python3

import glob
import os
import re
import sys
from argparse import Namespace

# https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit/tree/master/Colours
g_geodata_qualitative = [
    "#FF1F5B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E", "#F28522",
    "#A0B1BA", "#A6761D", "#E9002D", "#FFAA00", "#00B000"
]
g_colorscheme = g_geodata_qualitative

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except:
    pass

try:
    import scipy
    import scipy.interpolate
except ImportError:
    pass

try:
    import numpy as np
except ImportError:
    pass

kThin = False

def GetParams():
    import cycler
    params = {
        'font.size': 8,
        'figure.figsize': (3.2, 2.4),
        #'figure.autolayout': True,
        'axes.autolimit_mode': 'round_numbers',
        'axes.xmargin': 0,
        'axes.ymargin': 0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.prop_cycle': cycler.cycler(color=g_colorscheme),
        'lines.markersize': 3,
        'lines.linewidth': 1.5,
    }
    return params

def ApplyParams(plt):
    plt.rcParams.update(GetParams())


# natural sort
def natkey(s, _nsre=re.compile('([0-9]+)')):
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(_nsre, s)
    ]


def natsorted(v):
    return sorted(v, key=natkey)


def Glob(d, fld):
    return natsorted(glob.glob(os.path.join(d, "{:}*.dat".format(fld))))


# Writes message to stdout
# s: string
def Log(s, noeol=False):
    if not noeol:
        s += "\n"
    sys.stdout.write(s)
    sys.stdout.flush()


def HideAxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)


def PlotGrid(ax, x1n, y1n, c='0', alpha=0.2):
    ax.set_xticks(x1n)
    ax.set_yticks(y1n)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    ax.grid(True, lw=0.5, c=c, alpha=alpha)


def InitBasicFigure(field, grid=False):
    print("InitBasicFigure is deprecated. Use InitFigure")
    figsize = 3.2
    resx = 640
    dpi = resx / figsize
    fig = plt.figure(figsize=(resx / dpi, resx / dpi), dpi=dpi)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    lx, ly = 1., 1.
    ny, nx = field.shape
    hx, hy = lx / nx, ly / ny
    x1 = (0.5 + np.arange(nx)) * hx
    y1 = (0.5 + np.arange(ny)) * hy
    x1node = np.arange(nx + 1) * hx
    y1node = np.arange(ny + 1) * hy
    HideAxis(ax)
    if grid:
        PlotGrid(ax, x1node, y1node)
    meta = dict()
    meta['x1'] = x1
    meta['y1'] = y1
    meta['fig'] = fig
    return fig, ax, meta


def InitFigure(nx, ny, grid=False, resx=640, figsizex=3.2):
    dpi = resx / figsizex
    resy = resx * ny / nx
    lx, ly = 1., ny / nx
    fig = plt.figure(figsize=(resx / dpi, resy / dpi), dpi=dpi)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    ax.set_xlim(0, lx)
    ax.set_ylim(0, ly)
    hx, hy = lx / nx, ly / ny
    x1c = (0.5 + np.arange(nx)) * hx
    y1c = (0.5 + np.arange(ny)) * hy
    x1n = np.arange(nx + 1) * hx
    y1n = np.arange(ny + 1) * hy
    HideAxis(ax)
    if grid:
        PlotGrid(ax, x1n, y1n)
    meta = Namespace
    meta.lx = lx
    meta.ly = ly
    meta.hx = hx
    meta.hy = hy
    meta.x1c = x1c
    meta.y1c = y1c
    meta.x1n = x1n
    meta.y1n = y1n
    meta.fig = fig
    return fig, ax, meta


def GetSquareFigure(grid=False, field=None, aspect=1, colorbar=False):
    assert aspect == 1 or not colorbar, "conflicting options"
    if colorbar:
        aspect = 0.88
    figsize = 3.2
    resx = 640
    resy = 640 / aspect
    dpi = resx / figsize
    fig = plt.figure(figsize=(resx / dpi, resy / dpi), dpi=dpi)
    ax = plt.Axes(fig, [0., 0., 1.,  resy / resx])
    fig.add_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    HideAxis(ax)
    if grid:
        x1 = np.linspace(0, 1, field.shape[1] + 1, endpoint=True)
        y1 = np.linspace(0, 1, field.shape[0] + 1, endpoint=True)
        PlotGrid(ax, x1, y1)
    matplotlib.rcParams['svg.hashsalt'] = 123  # for reproducible svg
    return fig, ax

def AddColorBar(fig, vmin, vmax, cmap):
    cbar_ax = fig.add_axes([0.15, 0.07, 0.7, 0.05])
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    vv = [vmin, (vmin + vmax) * 0.5, vmax]
    if abs(vv[1]) < 1e-14:
        vv[1] = 0
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal', ticks=vv)
    cbar.ax.set_xticklabels(["{:.4g}".format(v) for v in vv])
    return cbar


def PlotFieldCoolwarm(ax, u, vmin=None, vmax=None):
    ax.imshow(np.flipud(u),
              vmin=vmin,
              vmax=vmax,
              extent=(0, 1, 0, 1),
              interpolation='nearest',
              cmap=plt.get_cmap("coolwarm"))


def PlotSquareField(ax, u, vmin=None, vmax=None, cmap="viridis"):
    ax.imshow(np.flipud(u),
              vmin=vmin,
              vmax=vmax,
              extent=(0, 1, 0, 1),
              interpolation='nearest',
              cmap=cmap)

def SaveFigure(fig, filename):
    matplotlib.rcParams['svg.hashsalt'] = 123  # for reproducible svg
    fig.savefig(filename)

def SaveBasicFigure(fig, filename):
    print("SaveBasicFigure is deprecated. Use SaveFigure")
    SaveFigure(fig, filename)



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
        if len(u.shape) > 2 and u.shape[2] == 2:
            u = u[:, :, :1]
        return u


# Converts field to image
def GetImg(u):
    s = u.shape
    assert len(s) in [2, 3]
    if len(s) == 3:
        u = u[:, :, 0]
    return np.flipud(u.T)


def PlotInitSq():
    fig, ax = plt.subplots(figsize=(3, 3))
    ax.set_aspect('equal')
    return fig, ax


def PlotInit():
    fig, ax = plt.subplots(figsize=(4, 3))
    return fig, ax


def PlotInit3():
    fig = plt.figure(figsize=(4, 3))
    ax = fig.gca(projection='3d')
    return fig, ax


def PlotSave(fig, ax, po):
    fig.tight_layout()
    fig.savefig(po, dpi=300)
    plt.close()


# u: 2d numpy array
# po: output path
def PlotFieldGray(ax, u, vmin=None, vmax=None):
    cmap = LinearSegmentedColormap.from_list('mycmap', ['#c8c5bd', '#787672'])
    ax.imshow(GetImg(u),
              extent=(0, 1, 0, 1),
              interpolation='nearest',
              vmin=vmin,
              vmax=vmax,
              cmap=cmap)


def PlotField(ax, field, meta, alpha=1, vmin=None, vmax=None, cmap=None):
    ax.imshow(np.flipud(field.T),
              extent=(0, meta.lx, 0, meta.ly),
              interpolation='nearest',
              alpha=alpha,
              vmin=vmin,
              vmax=vmax,
              cmap=None)

def PlotFieldText(ax, field, meta, skip_values=[], fmt="{:.2g}"):
    for ix, iy in np.ndindex(field.shape):
        value = field[ix, iy]
        if value in skip_values:
            continue
        ax.text(
            meta.x1c[ix],
            meta.y1c[iy],
            fmt.format(value),
            horizontalalignment='center',
            verticalalignment='center',
        )


# u: 2d numpy array
# po: output path
def PlotFieldBwr(ax, u, vmin=None, vmax=None):
    ax.imshow(GetImg(u),
              extent=(0, 1, 0, 1),
              interpolation='nearest',
              vmin=vmin,
              vmax=vmax,
              cmap=plt.get_cmap("bwr"))



def PlotLines(ax, xa, ya, xb, yb):
    xa = xa.flatten()
    ya = ya.flatten()
    xb = xb.flatten()
    yb = yb.flatten()
    xy = np.vstack((xa, xb, ya, yb))
    xy = xy.T
    xy = xy.reshape((xa.size * 2, 2))
    if kThin:
        ax.plot(*xy, c='k', alpha=0.9, lw=0.5)
    else:
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

    if len(e) == 1:  # if only one point found, set second to the same
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
    x, y, z, c = (d.T)[0:4, :]
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
            if kThin:
                ax.plot(x[ti], y[ti], c=cl, zorder=10, lw=0.5, alpha=0.5)
                ax.scatter(x[ti],
                           y[ti],
                           c=cl,
                           s=0.5,
                           lw=0.2,
                           zorder=11,
                           alpha=0.8,
                           edgecolor='black')
            else:
                ax.plot(x[ti], y[ti], c=cl, zorder=10, lw=1, alpha=0.5)
                ax.scatter(x[ti],
                           y[ti],
                           c=cl,
                           s=2,
                           lw=0.3,
                           zorder=11,
                           edgecolor='black')


def IsGerris(s):
    # check if gerris (odd nx)
    return s[0] % 2 == 1


# Returns dimension by shape
def GetDim(s):
    if len(s) == 2 or s[2] in [1, 2]:
        return 2
    return 3


# s: shape of data array
# Returns:
# x1,y1,z1: coordinates of cell centers
# hx,hy,hz: mesh steps (hz=1 if 2d))
def GetGeom(s):
    [nx, ny, nz] = s
    # cells
    ge = IsGerris(s)
    sx, sy, sz = (nx, ny, nz) if not ge else (nx - 1, ny - 1, nz - 1)
    ext = 1.  # extent
    h = ext / max(sx, sy, sz)
    hx = h
    hy = h
    hz = h
    # mesh
    if ge:
        x1 = (0. + np.arange(nx)) * hx
        y1 = (0. + np.arange(ny)) * hy
        z1 = (0. + np.arange(nz)) * hz
    else:
        x1 = (0.5 + np.arange(nx)) * hx
        y1 = (0.5 + np.arange(ny)) * hy
        z1 = (0.5 + np.arange(nz)) * hz
    return x1, y1, z1, hx, hy, hz


def GetMesh(x, y, z):
    return np.meshgrid(x, y, z, indexing='ij')


# assuming uniform mesh
def GetMeshStep(s):
    if IsGerris(s):
        return 1. / (s[0] - 1)
    return 1. / s[0]


# Returns mesh nodes
# assume [0,1]x[0,1]x[0,1] domain
def GetMeshNodes(hx, hy, hz):
    xn1 = np.arange(0., 1. + hx * 0.5, hx)
    yn1 = np.arange(0., 1. + hy * 0.5, hy)
    zn1 = np.arange(0., 1. + hz * 0.5, hz)
    return xn1, yn1, zn1


# Radius of circle form area
def GetReff2(a):
    return (a / np.pi)**0.5


# Radius of sphere from volume
def GetReff3(v):
    return (v / np.pi * 3. / 4.)**(1. / 3.)


# Returns effective radius
def GetReff(vf, dim):
    x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
    # cell volume
    vc = hx * hy * hz
    # integral of vf
    v = vf.sum() * vc
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


# Reads field by prefix
# pt: path template
# fld: field prefix
def ReadField(pt, fld):
    p = GetFieldPath(pt, fld)
    return ReadArray(GetFieldPath(pt, fld))


# Save figure with title only
# t: title
# po: output path
def FigTitle(t, po):
    fig, ax = plt.subplots(figsize=(5, 0.5))
    ax.text(0., 0., t, fontsize=15)
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    #fig.tight_layout()
    #fig.savefig(po, bbox_inches='tight')
    fig.savefig(po)
    plt.close()


# Figure with volume fraction
# pt: path template
# vfn: volume fraction field name
def FigVf(pt, vfn='vf'):
    vf = ReadField(pt, vfn)
    dim = GetDim(vf.shape)
    x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
    x, y, z = GetMesh(x1, y1, z1)
    ii = np.where(np.isfinite(vf))

    # bubble
    cx, cy, cz, rx, ry, rz = LoadBub()

    # center of mass
    #cx = np.sum(x[ii] * vf[ii]) / np.sum(vf[ii])
    #cy = np.sum(y[ii] * vf[ii]) / np.sum(vf[ii])
    #cz = np.sum(z[ii] * vf[ii]) / np.sum(vf[ii])

    # xy slice through center of mass
    iz = np.argmin(abs(cz - z))
    if vf is not None: vf = vf[:, :, iz:iz + 1]

    fig, ax = PlotInitSq()
    x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
    xn1, yn1, zn1 = GetMeshNodes(hx, hy, hz)
    x, y, z = GetMesh(x1, y1, z1)

    # grid
    PlotGrid(ax, xn1, yn1)
    # field
    PlotFieldGray(ax, vf, vmin=0., vmax=1.)
    # lines
    if dim == 2:
        a, nx, ny = [ReadField(pt, n) for n in ['a', 'nx', 'ny']]
        if all(u is not None for u in [a, nx, ny]):
            l = GetLines(x, y, a, nx, ny, hx, hy, vf)
            PlotLines(ax, *l)
    # partilces
    pa = GetFieldPath(pt, "partit")
    if os.path.isfile(pa) and dim == 2:
        Log(pa)
        PlotPart(ax, pa, sk=1)
    # closeup
    r = max(rx, ry) * 1.2
    ax.set_xlim([cx - r, cx + r])
    ax.set_ylim([cy - r, cy + r])
    # save
    po = GetFieldPath(pt, "vf", "pdf")
    Log(po)
    PlotSave(fig, ax, po)


# Exact trajectory.
# x0: initial center, shape (3)
# vx0: velocity, shape (3)
# tmax: total time
# n: number of segments
def GetTrajE(x0, vx0, tmax, n=150):
    x = [np.linspace(x0[d], x0[d] + vx0[d] * tmax, n) for d in range(3)]
    return x


# Extracts trajectory of center of mass.
# pp: list of paths to volume fraction fields
# Returns:
# x,y: arrays
def GetTraj(pp):
    # result
    cx = []
    cy = []
    cz = []
    for i, p in enumerate(pp):
        # report
        Log("{:} ".format(i), True)
        # read array
        vf = ReadArray(p)
        x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
        x, y, z = GetMesh(x1, y1, z1)
        # compute center
        cx0 = (x * vf).sum() / vf.sum()
        cy0 = (y * vf).sum() / vf.sum()
        cz0 = (z * vf).sum() / vf.sum()
        cx.append(cx0)
        cy.append(cy0)
        cz.append(cz0)
    Log("")

    cx = np.array(cx)
    cy = np.array(cy)
    cz = np.array(cz)
    return cx, cy, cz


# Extracts value of field along trajectory.
# pp: list of paths to volume fraction fields
# cx,cy,cz: trajectory
# fld: field prefix (e.g. 'p')
# Returns:
# cu: field at center
def GetTrajFld(pp, cx, cy, cz, fld):
    cu = []
    for i, p in enumerate(pp):
        # report
        Log("{:} ".format(i), True)
        # path template
        pt = GetPathTemplate(p)
        # read field
        u = ReadField(pt, fld)
        # mesh
        x1, y1, z1, hx, hy, hz = GetGeom(u.shape)
        x, y, z = GetMesh(x1, y1, z1)
        # index
        ci0 = int(cx[i] / hx)
        cj0 = int(cy[i] / hy)
        ck0 = int(cz[i] / hz)
        # field
        cu.append(u[ci0, cj0, ck0])
    Log("")

    cu = np.array(cu)
    return cu


# Average field weighted with volume fraction.
# p: path to volume fraction
# fld: field prefix (e.g. 'p')
# Returns:
# cu: average value of field
def GetAvgFld0(p, fld):
    # read volume fraction
    v = ReadArray(p)
    # read field
    u = ReadField(GetPathTemplate(p), fld)
    # field value
    return (u * v).mean() / v.mean()


# pp: list of path to volume fraction
# fld: field prefix
def GetAvgFld(pp, fld):
    return np.array([GetAvgFld0(p, fld) for p in pp])


# Function applied to field.
# f: function
# p: path to volume fraction
# fld: field prefix (e.g. 'p')
# Returns:
# cu: f(u)
def GetFuncFld0(f, p, fld):
    u = ReadField(GetPathTemplate(p), fld)
    return f(u)


# pp: list of path to volume fraction
# fld: field prefix
def GetFuncFld(f, pp, fld):
    return np.array([GetFuncFld0(f, p, fld) for p in pp])


# Magnitude of difference from uniform reference field.
# pp: list of paths to volume fraction fields
# fld: field prefix (e.g. 'vx')
# ue: reference
# Returns:
# umax: max norm
# u1: l1 norm
# u2: l2 norm
def GetDiff(pp, fld, ue):
    umax = []
    u1 = []
    u2 = []
    for i, p in enumerate(pp):
        # read volume fraction
        v = ReadArray(p)
        # read field
        u = ReadField(GetPathTemplate(p), fld)
        # difference
        du = abs(u - ue)
        # norms
        umax.append((du * v).max())
        u1.append((du * v).sum() / v.sum())
        u2.append(((du**2 * v).sum() / v.sum())**0.5)
    umax = np.array(umax)
    u1 = np.array(u1)
    u2 = np.array(u2)
    return umax, u1, u2


# Writes line as x,y columns.
# x,y: 1d arrays
# po: output path
def WriteLine0(x, y, po):
    np.savetxt(po, np.array((x, y)).T)


# Output path for single line.
# l: line label
# po: output path [base].[ext]
# Returns:
# "[base]_[l].dat"
def GetLinePath(l, po):
    b = os.path.splitext(po)[0]
    return "{:}_{:}.dat".format(b, l)


# Writes line as x,y columns.
# x,y: 1d arrays
# l: line label
# po: output path [base].[ext]
# Output:
# [base]_[l].dat: file with x,y columns
def WriteLine(x, y, l, po):
    WriteLine0(x, y, GetLinePath(l, po))


# Plots trajectories
# xx,yy: list of arrays for axes
# ll: line labels
# lx,ly: axes labels
# ylog: log-scale in y
# po: output path
def PlotTrajFld(xx,
                yy,
                ll,
                lx,
                ly,
                po,
                vmin=None,
                vmax=None,
                ystep=None,
                ylog=False):
    fig, ax = PlotInit()
    if vmin is None: vmin = yy[0].min()
    if vmax is None: vmax = yy[0].max()
    i = 0
    for x, y, l in zip(xx, yy, ll):
        WriteLine(x, y, l, po)
        if i == len(xx) - 1:  # separate for last
            ax.plot(x, y, label=l, c="0.5", ls='--')
        else:
            ax.plot(x, y, label=l)
        i += 1
    ax.legend()
    ax.grid(True)
    ax.set_xlabel(lx)
    ax.set_ylabel(ly)
    ax.set_ylim(vmin, vmax)
    if ylog:
        ax.set_yscale('log')
    if ystep is not None:
        ax.set_yticks(np.arange(vmin, vmax + 1e-10, ystep))
    PlotSave(fig, ax, po)


# Plots trajectories
# xx,yy: list of arrays
# ll: labels
# s: shape of of data array
# po: output path
def PlotTraj(xx, yy, ll, s, po):
    fig, ax = PlotInitSq()
    ax.set_aspect('equal')
    i = 0
    for x, y, l in zip(xx, yy, ll):
        WriteLine(x, y, l, po)
        if i == len(xx) - 1:  # different style for last (exact)
            ax.plot(x, y, label=l, c="0.5", ls='--')
        else:
            ax.plot(x, y, label=l)
        i += 1
    ax.legend()
    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    x1, y1, z1, hx, hy, hz = GetGeom(s)
    xn1, yn1, zn1 = GetMeshNodes(hx, hy, hz)
    PlotGrid(ax, xn1, yn1)

    q, q, q, rx, ry, rz = LoadBub()
    x1, y1, z1, hx, hy, hz = GetGeom(s)
    plt.title("r/h={:0.3f}".format(rx / hx))
    PlotSave(fig, ax, po)


# Computes norms: max, L1, L2
# u: array
# Returns:
# m, l1, l2
def GetNorm(k):
    k = abs(k)
    m = k.max()
    l1 = k.mean()
    l2 = (k**2).mean()**0.5
    return m, l1, l2


# bubble location
def LoadBub():
    bb = np.loadtxt("b.dat")
    if len(bb.shape) == 2:
        return bb[0, :]
    return bb


def GetConf(p):
    c = {}
    exec(open(p).read(), None, c)
    return c


# initial velocity
def LoadVel():
    return np.loadtxt("vel")


# surface tension
def LoadSig():
    return np.loadtxt("sig")


# viscosity
def LoadMu():
    return np.loadtxt("mu")


# viscosity
def LoadTmax():
    return np.loadtxt("tmax")


# Returns mean curvature of surface z=h(x,y) at point (x,y,z)
# h: function h(x,y)
# x,y: coordinates
# dx: step to compute derivative
def GetCurv(h, x, y, dx=1e-3):
    dy = dx

    # returns df/dx of function f(x,y)
    def Dx(f):
        return lambda x, y: (f(x + dx, y) - f(x - dx, y)) / (2. * dx)

    # returns df/dy of function f(x,y)
    def Dy(f):
        return lambda x, y: (f(x, y + dy) - f(x, y - dy)) / (2. * dy)

    hx = Dx(h)(x, y)
    hy = Dy(h)(x, y)
    hxx = Dx(Dx(h))(x, y)
    hxy = 0.5 * (Dx(Dy(h))(x, y) + Dy(Dx(h))(x, y))
    hyy = Dy(Dy(h))(x, y)

    #a = (1. + hx ** 2) * hyy - 2. * hx * hy * hxy + (1. + hy ** 2) * hxx
    #c = (hx ** 2 + hy ** 2 + 1.) ** (3. / 2)
    #return -a / c

    a = hx**2 * hxx + 2. * hx * hy * hxy + hy**2 * hyy
    b = (hx**2 + hy**2 + 1.) * (hxx + hyy)
    c = (hx**2 + hy**2 + 1.)**(3. / 2)
    return (a - b) / c


# Prints stat of field x with label lbl.
def P(x, lbl):
    print("{:}: shape={:} min={:}, max={:}, avg={:}".format(
        lbl, x.shape, x.min(), x.max(), x.mean()))


# Returns curvature of ellipse
# (x/rx)**2 + (y/ry)**2 = 1 at point x
def GetEllip2Curv0(x, rx, ry):
    if rx == ry:
        return 1. / rx
    dx = rx * 1e-3
    x = np.clip(x, -rx + dx * 2, rx - dx * 2)  # XXX clip
    return GetCurv(lambda xx, zz: max(0., 1. - (xx / rx)**2)**0.5 * ry,
                   x,
                   0.,
                   dx=dx)


# Returns curvature of ellipse
# (x/rx)**2 + (y/ry)**2 = 1 at point x,y
def GetEllip2Curv(x, y, rx, ry):
    if abs(x) < abs(y):
        return GetEllip2Curv0(x, rx, ry)
    return GetEllip2Curv0(y, ry, rx)


# Returns curvature of ellipsoid
# (x/rx)**2 + (y/ry)**2 + (z/rz) ** 2 = 1 at point x,y
def GetEllip3Curv0(x, y, rx, ry, rz):
    if rx == ry and ry == rz:
        return 2. / rx
    dd = 1e-3
    dx = rx * dd
    dy = ry * dd
    return GetCurv(lambda xx, yy: max(0., 1. - (xx / rx)**2 -
                                      (yy / ry)**2)**0.5 * rz,
                   x,
                   y,
                   dx=dx)


# Returns curvature of ellipsoid
# (x/rx)**2 + (y/ry)**2 + (z/rz) ** 2 = 1 at point x,y,z
def GetEllip3Curv(x, y, z, rx, ry, rz):
    if abs(z) >= abs(x) and abs(z) >= abs(y):  # z(x,y)
        return GetEllip3Curv0(x, y, rx, ry, rz)
    elif abs(y) >= abs(x) and abs(y) >= abs(z):  # y(x,z)
        return GetEllip3Curv0(x, z, rx, rz, ry)
    # x(y,z)
    return GetEllip3Curv0(y, z, ry, rz, rx)


# Evaluates exact curvature of ellipsoid
# dim: dimension, 2 or 3
# x,y,z: points
# Returns:
# k: curvature of shape x
def GetExactK(dim, x, y, z):
    cx, cy, cz, rx, ry, rz = LoadBub()

    dx = x - cx
    dy = y - cy
    dz = z - cz

    k = np.zeros_like(x)

    if dim == 2:
        u = np.arctan2(dy, dx)
        r = ((np.cos(u) / rx)**2 + (np.sin(u) / ry)**2)**(-0.5)
        dx = r * np.cos(u)
        dy = r * np.sin(u)
        for i, wx, wy in zip(range(len(dx)), dx, dy):
            k[i] = GetEllip2Curv(wx, wy, rx, ry)

    if dim == 3:
        u = np.arctan2(dy, dx)
        v = np.arctan2(dz, (dx**2 + dy**2)**0.5)
        r = ((np.cos(u) * np.cos(v) / rx)**2 +
             (np.sin(u) * np.cos(v) / ry)**2 + (np.sin(v) / rz)**2)**(-0.5)
        dx = r * np.cos(u) * np.cos(v)
        dy = r * np.sin(u) * np.cos(v)
        dz = r * np.sin(v)
        for i, wx, wy, wz in zip(range(len(dx)), dx, dy, dz):
            k[i] = GetEllip3Curv(wx, wy, wz, rx, ry, rz)

    return k


# Histogram of curvature
def FigHistK(vf, kk, ll, po, title=None):
    dim = GetDim(vf.shape)
    x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
    x, y, z = GetMesh(x1, y1, z1)

    # interface cells
    th = 1e-8
    ii = np.where((vf > th) & (vf < 1. - th))
    x = x[ii]
    y = y[ii]
    z = z[ii]

    # exact curvature
    ke = GetExactK(dim, x, y, z)
    # average curvature
    kea = ke.mean()

    def ER(k):
        if IsGerris(vf.shape):
            k = -k
        k = k[ii]
        return (k - ke) / kea  # error

    fig, ax = PlotInit()
    for k, l in zip(kk, ll):
        er = ER(k)
        rm = 1.
        h, b = np.histogram(er, 200, range=(-rm, rm), density=False)
        bc = (b[1:] + b[:-1]) * 0.5
        WriteLine(bc, h, l, po)
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
    cx, cy, cz, rx, ry, rz = LoadBub()
    with open(GetLinePath("norm", po), 'w') as f:
        for k, l in zip(kk, ll):
            if k is not None:
                m, l1, l2 = GetNorm(ER(k))
                f.write("{:} {:} {:} {:} {:} {:}\n".format(
                    l, rx / hx, ry / hx, m, l1, l2))


# Figure of curvature vs angle.
# vf: volume fraction field
# kk: list of curvature fields of shape vf.shape
# ll: list of labels
def FigAng(vf, kk, ll, po):
    dim = GetDim(vf.shape)
    x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
    x, y, z = GetMesh(x1, y1, z1)

    # select interface points from x,y,z
    th = 0.
    ii = np.where((vf > th) & (vf < 1. - th))
    x = x[ii]
    y = y[ii]
    z = z[ii]

    cx, cy, cz, rx, ry, rz = LoadBub()

    # cartesian to spherical for points x,y,z on interface
    dx = x - cx
    dy = y - cy
    dz = z - cz
    u = np.arctan2(dy, dx)
    v = np.arctan2(dz, (dx**2 + dy**2)**0.5)
    # angles of interface points in deg
    degu = np.degrees(u)
    degv = np.degrees(v)

    # hires angle
    uh = np.linspace(-np.pi, np.pi, 200)
    if dim == 3:
        vh = np.linspace(-np.pi * 0.5, np.pi * 0.5, 200)
    else:
        vh = np.array([0.])
    # hires points
    # ellipsoid in spherical coordinates: r=r(u,v)
    r = ((np.cos(uh) * np.cos(vh) / rx)**2 +
         (np.sin(uh) * np.cos(vh) / ry)**2 + (np.sin(vh) / rz)**2)**(-0.5)
    xh = cx + r * np.cos(uh) * np.cos(vh)
    yh = cy + r * np.sin(uh) * np.cos(vh)
    zh = cz + r * np.sin(vh)
    # hires curvature
    keh = GetExactK(dim, xh.flatten(), yh.flatten(), zh.flatten())
    keh = keh.reshape(xh.shape)
    # average curvature
    kea = keh.mean()

    # line in angles
    nls = 200
    us = np.linspace(-np.pi, np.pi, nls)
    if dim == 3:
        vs = np.linspace(-np.pi * 0.5, np.pi * 0.5, nls)
    else:
        vs = np.linspace(0., 0., nls)
    # angles of line in deg
    degus = np.degrees(us)
    degvs = np.degrees(vs)

    # Interpolate to line
    # u,v: source points
    # k: source field
    # ut,vt: target points
    # peru: period in u
    # Returns:
    # kt: field k interpolated from u,x to ut,vt
    def IL(u, v, k, ut, vt, peru=None):
        if v.ptp() > 0.:  # 3d
            return scipy.interpolate.griddata((u, v),
                                              k, (ut, vt),
                                              method='nearest')
        else:  # 2d
            if peru is not None:
                u = np.hstack([u - peru, u, u + peru])
                k = np.hstack([k, k, k])
            fk = scipy.interpolate.interp1d(u,
                                            k,
                                            bounds_error=False,
                                            kind='nearest')
            return fk(ut)

    # interpolate curvature to line
    kks = []
    for k in kk:
        if k is not None:
            kks.append(IL(degu, degv, k[ii], degus, degvs, peru=360.))
        else:
            kks.append(None)

    fig, ax = PlotInit()

    # plot interpolated curvature along line
    for k, l in zip(kks, ll):
        if k is not None:
            WriteLine(degus, k / kea, l, po)
            ax.plot(degus, k / kea, label=l)

    # plot lowres exact curvature
    ke = GetExactK(dim, x, y, z)
    kes = IL(degu, degv, ke, degus, degvs, peru=360.)
    WriteLine(degus, kes / kea, "exact", po)
    ax.plot(degus, kes / kea, label="exact")

    # plot hires exact curvature
    WriteLine(degus, keh / kea, "exacth", po)
    ax.plot(degus, keh / kea, label="exact", c="0.5", ls='--')

    ax.set_xlabel(r"$\varphi$ [deg]")
    ax.set_ylabel(r"normalized curvature")
    ax.grid()
    ax.legend()
    ax.set_ylim(0., 2.)
    PlotSave(fig, ax, po)

    poe = 'er' + po

    fig, ax = PlotInit()
    # plot interpolated curvature along line
    ke = GetExactK(dim, x, y, z)
    kes = IL(degu, degv, ke, degus, degvs, peru=360.)
    for k, l in zip(kks, ll):
        if k is not None:
            WriteLine(degus, (k - kes) / kea, l, poe)
            ax.plot(degus, abs((k - kes) / kea), label=l)

    ax.set_xlabel(r"$\varphi$ [deg]")
    ax.set_ylabel(r"curvature error")
    ax.set_yscale('log')
    ax.set_ylim(1e-5, 3)
    ax.grid()
    ax.legend()
    PlotSave(fig, ax, poe)


# Figure with curvature depending on angle, 2D
# XXX: deprecated
# vf: volume fraction field
# kk: list of curvature fields
# ll: list of labels
def FigAng2(vf, kk, ll, po):
    dim = GetDim(vf.shape)
    x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
    x, y, z = GetMesh(x1, y1, z1)

    # interface cells
    th = 0
    ii = np.where((vf > th) & (vf < 1. - th))
    x = x[ii]
    y = y[ii]
    z = z[ii]

    # hires angle
    anh = np.linspace(-np.pi, np.pi, 200)
    degh = np.degrees(anh)
    # hires points
    cx, cy, cz, rx, ry, rz = LoadBub()
    r = (rx * ry) / ((ry * np.cos(anh))**2 + (rx * np.sin(anh))**2)**0.5
    xh = cx + r * np.cos(anh)
    yh = cy + r * np.sin(anh)
    # hires curvature
    keh = GetExactK(dim, xh, yh, cz)
    # average curvature
    kea = keh.mean()

    dx = x - cx
    dy = y - cy
    an = np.arctan2(dy, dx)
    deg = np.degrees(an)
    s = np.argsort(an)

    fig, ax = PlotInit()

    # plot estimates
    for k, l in zip(kk, ll):
        if k is not None:
            ax.plot(deg[s], k[ii][s] / kea, label=l)

    # plot lowres exact curvature
    ax.plot(deg[s], GetExactK(dim, x, y, z)[s] / kea, label="exact")

    # plot exact curvature
    ax.plot(degh, keh / kea, label="exact", c="0.5", ls='--')

    ax.set_xlabel(r"angle [deg]")
    ax.set_ylabel(r"normalized curvature")
    q = 180.
    ax.set_xticks(np.arange(-q, q * 1.1, 90.))
    ax.set_xlim(-q, q)
    ax.legend()
    ax.set_ylim(0., 2.)
    ax.grid()
    PlotSave(fig, ax, po)


def Univel():
    # directories
    dd = ['ch', 'ge', 'ba']
    # labels
    ll = []
    # trajectories
    xx = []
    yy = []
    zz = []
    # pressure along trajectories
    ee = []
    # velocity error max
    vvxm = []
    vvym = []
    vvzm = []
    # velocity error l2
    vvx2 = []
    vvy2 = []
    vvz2 = []
    # average velocity
    vvx = []
    vvy = []
    vvz = []

    # volume fraction from dd[0]
    pp = Glob(dd[0], "vf")
    pt = GetPathTemplate(pp[0])
    vf = ReadField(pt, "vf")

    dim = GetDim(vf.shape)

    conf = GetConf("./par.py")

    cx, cy, cz, rx, ry, rz = LoadBub()
    x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
    sig = conf['sig']
    mu = conf['mu']

    # exact velocity
    vele = np.array(conf['vel'])
    # exact trajectories
    tmax = conf['tmax']
    x, y, z = GetTrajE([cx, cy, cz], vele, tmax)
    zs = x * 0.  # zeros
    # exact average curvature
    if dim == 3:
        vol = 4. / 3. * np.pi * rx * ry * rz
        req = (vol * 3. / (4. * np.pi))**(1. / 3.)
        kavg = 2. / req
    else:
        vol = np.pi * rx * ry
        req = (vol / np.pi)**0.5
        kavg = 1. / req
    # exact pressure jump
    eex = sig * kavg
    e = zs + 1
    # error in velocity
    vx = x * 0
    vy = x * 0
    vz = x * 0
    # append
    xx.append(x)
    yy.append(y)
    zz.append(z)
    vvxm.append(vx)
    vvym.append(vy)
    vvzm.append(vz)
    vvx2.append(vx)
    vvy2.append(vy)
    vvz2.append(vz)
    vvx.append(vx)
    vvy.append(vy)
    vvz.append(vz)
    ee.append(e)
    ll.append("exact")

    # lines
    ll.append('ch')
    ll.append('ge')
    ll.append('ba')
    for d, l in zip(dd, ll):
        print(d)
        pp = Glob(d, "vf")
        if len(pp) == 0:
            continue
        pt = GetPathTemplate(pp[0])
        # volume fraction
        FigVf(pt)
        # trajectory
        x, y, z = GetTraj(pp)
        # pressure jump relative to exact
        e = GetFuncFld(np.ptp, pp, 'p') / eex
        # error in velocity
        vxm, vx1, vx2 = GetDiff(pp, "vx", vele[0])
        vym, vy1, vy2 = GetDiff(pp, "vy", vele[1])
        # velocity
        vx = GetAvgFld(pp, "vx")
        vy = GetAvgFld(pp, "vy")
        # velocity dimensionless factor
        q = mu / sig
        q = 1.  # XXX
        if dim == 3:
            vzm, vz1, vz2 = GetDiff(pp, "vz", vele[2])
            vz = GetAvgFld(pp, "vz")
            vvzm.append(vzm * q)
            vvz2.append(vz2 * q)
            vvz.append(vz)
        xx.append(x)
        yy.append(y)
        zz.append(z)
        ee.append(e)
        vvxm.append(vxm * q)
        vvym.append(vym * q)
        vvx2.append(vx2 * q)
        vvy2.append(vy2 * q)
        vvx.append(vx)
        vvy.append(vy)

    # reorder lines
    xx = xx[1:] + [xx[0]]
    yy = yy[1:] + [yy[0]]
    zz = zz[1:] + [zz[0]]
    ee = ee[1:] + [ee[0]]
    vvxm = vvxm[1:] + [vvxm[0]]
    vvym = vvym[1:] + [vvym[0]]
    vvx2 = vvx2[1:] + [vvx2[0]]
    vvy2 = vvy2[1:] + [vvy2[0]]
    vvx = vvx[1:] + [vvx[0]]
    vvy = vvy[1:] + [vvy[0]]
    if dim == 3:
        vvzm = vvzm[1:] + [vvzm[0]]
        vvz2 = vvz2[1:] + [vvz2[0]]
        vvz = vvz[1:] + [vvz[0]]
    ll = ll[1:] + [ll[0]]

    # time
    tt = []
    for i, x in enumerate(xx):
        tt.append(np.linspace(0., conf['tmax'], len(x)))

    # Plot trajectories
    PlotTraj(xx, yy, ll, vf.shape, "traj.pdf")
    # Plot pressure jump
    PlotTrajFld(tt, ee, ll, "t", "p", "trajp.pdf", vmin=0., vmax=2.)
    # Plot x,y
    gp = 0.05  # gap
    PlotTrajFld(tt,
                xx,
                ll,
                "t",
                "x",
                "trajx.pdf",
                ystep=0.1,
                vmin=cx - gp,
                vmax=cx + vele[0] * tmax + gp)
    PlotTrajFld(tt,
                yy,
                ll,
                "t",
                "y",
                "trajy.pdf",
                ystep=0.1,
                vmin=cy - gp,
                vmax=cy + vele[1] * tmax + gp)
    if dim == 3:
        PlotTrajFld(tt,
                    zz,
                    ll,
                    "t",
                    "z",
                    "trajz.pdf",
                    ystep=0.1,
                    vmin=cz - gp,
                    vmax=cz + vele[2] * tmax + gp)
    # Plot error in velocity
    vmin = 1e-12
    vmax = 1e1
    PlotTrajFld(tt,
                vvxm,
                ll,
                "t",
                "vx",
                "trajvxm.pdf",
                vmin=vmin,
                vmax=vmax,
                ylog=True)
    PlotTrajFld(tt,
                vvym,
                ll,
                "t",
                "vy",
                "trajvym.pdf",
                vmin=vmin,
                vmax=vmax,
                ylog=True)
    PlotTrajFld(tt,
                vvx2,
                ll,
                "t",
                "vx",
                "trajvx2.pdf",
                vmin=vmin,
                vmax=vmax,
                ylog=True)
    PlotTrajFld(tt,
                vvy2,
                ll,
                "t",
                "vy",
                "trajvy2.pdf",
                vmin=vmin,
                vmax=vmax,
                ylog=True)
    PlotTrajFld(tt, vvx, ll, "t", "vx", "trajvx.pdf")
    PlotTrajFld(tt, vvy, ll, "t", "vy", "trajvy.pdf")
    if dim == 3:
        PlotTrajFld(tt,
                    vvzm,
                    ll,
                    "t",
                    "vz",
                    "trajvzm.pdf",
                    vmin=vmin,
                    vmax=vmax,
                    ylog=True)
        PlotTrajFld(tt,
                    vvz2,
                    ll,
                    "t",
                    "vz",
                    "trajvz2.pdf",
                    vmin=vmin,
                    vmax=vmax,
                    ylog=True)
        PlotTrajFld(tt, vvz, ll, "t", "vz", "trajvz.pdf")


def Curv():
    # directories
    dd = ['ch', 'ge']
    # labels
    ll = ['mfer', 'gerris']
    # curvature arrays
    kk = []

    for d, l in zip(dd, ll):
        print(d)
        pp = Glob(d, 'u')
        assert pp is not None
        pt = GetPathTemplate(pp[-1])

        # volume fraction
        FigVf(pt, 'u')

        # append curvature
        k = ReadField(pt, 'k')
        assert k is not None
        if IsGerris(k.shape):
            k = k * (-1)
        kk.append(k)

    # volume fraction from dd[0]
    pp = Glob(dd[0], 'u')
    pt = GetPathTemplate(pp[-1])
    vf = ReadField(pt, "u")
    dim = GetDim(vf.shape)

    # title
    cx, cy, cz, rx, ry, rz = LoadBub()
    x1, y1, z1, hx, hy, hz = GetGeom(vf.shape)
    if dim == 2:
        title = "rx/h={:0.2f} ry/h={:0.2f}".format(rx / hx, ry / hx)
    else:
        title = "rx/h={:0.2f} ry/h={:0.2f} rz/h={:0.2f}".format(
            rx / hx, ry / hx, rz / hx)

    # curvature vs angle
    po = 'ang.pdf'
    FigAng(vf, kk, ll, po)

    # curvature histogram
    po = 'hist.pdf'
    FigHistK(vf, kk, ll, po, title=title)

def GetStep(path):
    return re.findall('[^_]*_([0-9]*)\.*.', os.path.basename(path))[0]

def GetSteps(paths):
    return list(map(GetStep, paths))

def ReplaceFilename(paths, pattern, keep_dir=True):
    """
    Replaces filename by pattern with step index substitution.
    paths: `list(str)`
        Paths.
    pattern: `str`
        Pattern containing a single `{}` to be replaced by step index.

    Example:
    >>> ReplaceFilename(["dir/vf_0001.xmf"], "sm_{}.vtk")
    >>> ["dir/sm_0001.vtk"]
    """
    r = []
    for f in paths:
        dirname = os.path.dirname(f)
        basename = os.path.basename(f)
        step = GetStep(f)
        if keep_dir:
            r.append(os.path.join(dirname, pattern.format(step)))
        else:
            r.append(pattern.format(step))
    return r
