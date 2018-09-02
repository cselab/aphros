#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re
import sys
from plotlib import *

def Log(s, noeol=False):
    if not noeol:
        s += s + "\n"
    sys.stdout.write(s)
    sys.stdout.flush()

# Figure with volume fraction
# pt: path template
# suf: output name suffix
def FigVf(pt):
    vf = ReadField(pt, 'vf')

    fig, ax = PlotInitSq()
    # q: dummy
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    xn1,yn1,zn1 = GetMeshNodes(hx, hy, hz)
    x,y,z = GetMesh(x1, y1, z1)

    # grid
    PlotGrid(ax, xn1, yn1)
    # field
    PlotFieldGray(ax, vf, vmin=0., vmax=1.)
    # partilces
    pa = GetFieldPath(pt, "partit")
    if os.path.isfile(pa):
        Log(pa)
        PlotPart(ax, pa, sk=4)

    # save
    po = GetFieldPath(pt, "vf", "pdf")
    Log(po)
    PlotSave(fig, ax, po)


# Exact trajectory.
# x0,y0: initial center
# vx0,vy0: velocity
# t: time
# n: frames
def GetTrajE(x0, y0, vx0, vy0, t, n=150):
    x = np.linspace(x0, x0 + vx0 * t, n)
    y = np.linspace(y0, y0 + vy0 * t, n)
    return x,y

# Extracts trajectory of center of mass.
# pp: list of paths to volume fraction fields
# Returns:
# x,y: arrays
def GetTraj(pp):
    # result
    cx = [] ; cy = []
    for i,p in enumerate(pp):
        # read array
        vf = ReadArray(p)
        x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
        x,y,z = GetMesh(x1, y1, z1)
        # compute center
        cx0 = (x * vf).sum() / vf.sum()
        cy0 = (y * vf).sum() / vf.sum()
        cx.append(cx0)
        cy.append(cy0)

    cx = np.array(cx)
    cy = np.array(cy)
    return cx,cy

# Average field weighed with volume fraction along trajectory.
# pp: list of paths to volume fraction fields
# fld: field prefix (e.g. 'p')
# Returns:
# cu: average field at every point in pp
def GetAvgFld(pp, fld):
    cu = []
    for i,p in enumerate(pp):
        # read field
        u = ReadField(GetPathTemplate(p), fld)
        # read volume fraction
        v = ReadArray(p)
        # field value
        cu.append((u * v).mean() / v.mean())
    cu = np.array(cu)
    return cu

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
    for i,p in enumerate(pp):
        # read field
        u = ReadField(GetPathTemplate(p), fld)
        # difference
        du = abs(u - ue)
        # norms
        umax.append(du.max())
        u1.append(du.mean())
        u2.append((du ** 2).mean() ** 0.5)
    umax = np.array(umax)
    u1 = np.array(u1)
    u2 = np.array(u2)
    return umax, u1, u2

# Extracts value of field along trajectory.
# pp: list of paths to volume fraction fields
# cx,cy: trajectory
# fld: field prefix (e.g. 'p')
# Returns:
# cu: field at center
def GetTrajFld(pp, cx, cy, fld):
    cu = []
    for i,p in enumerate(pp):
        # read field
        u = ReadField(GetPathTemplate(p), fld)
        # mesh
        x1,y1,z1,hx,hy,hz = GetGeom(u.shape)
        x,y,z = GetMesh(x1, y1, z1)
        # index
        ci0 = int(cx[i] / hx)
        cj0 = int(cy[i] / hy)
        # field value
        cu.append(u[ci0, cj0, 0])
    cu = np.array(cu)
    return cu

# Plots field along trajectory
# xx,yy: list of arrays for axes
# ll: line labels
# lx,ly: axes labels
# ylog: log-scale in y
# po: output path
def PlotTrajFld(xx, yy, ll, lx, ly, po, vmin=None, vmax=None, ystep=None,
                ylog=False):
    global reff, hx
    fig, ax = PlotInit()
    if vmin is None: vmin = yy[0].min()
    if vmax is None: vmax = yy[0].max()
    i = 0
    for x,y,l in zip(xx, yy, ll):
        if i == len(xx) - 1: # custom for last
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
# xx,yy: list of arrays with coordinates
# ll: labels
# s: shape of of data array
def PlotTraj(xx, yy, ll, s):
    global reff, hx
    fig, ax = PlotInitSq()
    ax.set_aspect('equal')
    i = 0
    for x,y,l in zip(xx, yy, ll):
        if i == len(xx) - 1: # custom for last
            ax.plot(x, y, label=l, c="0.5", ls='--')
        else:
            ax.plot(x, y, label=l)
        i += 1
    ax.legend()
    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    x1,y1,z1,hx,hy,hz = GetGeom(s)
    xn1,yn1,zn1 = GetMeshNodes(hx, hy, hz)
    PlotGrid(ax, xn1, yn1)

    q,q,q,rx,ry,rz = LoadBub()
    x1,y1,z1,hx,hy,hz = GetGeom(s)
    plt.title("r/h={:0.3f}".format(rx / hx))
    po = 'traj.pdf'
    PlotSave(fig, ax, po)

def Main():
    # directories
    dd = ['ch', 'ge']
    # labels
    ll = []
    # trajectories
    xx = []
    yy = []
    # pressure jump error
    ee = []
    # velocity error max
    vvxm = []
    vvym = []
    vvzm = []
    # velocity error l2
    vvx2 = []
    vvy2 = []
    vvz2 = []

    # volume fraction from dd[0]
    pp = Glob(dd[0], "vf")
    pt = GetPathTemplate(pp[0])
    vf = ReadField(pt, "vf")

    cx,cy,cz,rx,ry,rz = LoadBub()
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)

    #vele = [0.4, 0.3, 0.]
    vele = [0., 0., 0.] # XXX

    # exact trajectories
    x,y = GetTrajE(cx, cy, vele[0], vele[1], 1.)
    # exact pressure jump
    sig = 1.
    eex = sig / rx
    # error in pressure jump
    e = (x * 0 + eex) / eex
    # error in velocity
    vx = x * 0
    vy = x * 0
    vz = x * 0
    # append
    xx.append(x)
    yy.append(y)
    ee.append(e)
    vvxm.append(vx)
    vvym.append(vy)
    vvzm.append(vz)
    vvx2.append(vx)
    vvy2.append(vy)
    vvz2.append(vz)
    ll.append("exact")

    # lines
    ll.append('mfer')
    ll.append('gerris')
    for d,l in zip(dd, ll):
        Log(d)
        pp = Glob(d, "vf")
        assert len(pp)
        pt = GetPathTemplate(pp[0])
        # volume fraction
        FigVf(pt)
        # trajectory
        x,y = GetTraj(pp)
        # pressure
        e = GetTrajFld(pp, x, y, "p")
        # pressure at corner
        e0 = GetTrajFld(pp, x * 0, y * 0, "p")
        # pressure jump relative to exact
        e = (e - e0) / eex
        # error in velocity
        vxm, vx1, vx2 = GetDiff(pp, "vx", vele[0])
        vym, vy1, vy2 = GetDiff(pp, "vy", vele[1])
        #vzmax, vz1, vz2 = GetDiff(pp, "vz", vele[2])
        # append
        xx.append(x)
        yy.append(y)
        ee.append(e)
        vvxm.append(vxm)
        vvym.append(vym)
        vvzm.append(vym) # TODO: vzmax
        vvx2.append(vx2)
        vvy2.append(vy2)
        vvz2.append(vy2) # TODO: vz2

    # reorder lines
    xx = xx[1:] + [xx[0]]
    yy = yy[1:] + [yy[0]]
    ee = ee[1:] + [ee[0]]
    vvxm = vvxm[1:] + [vvxm[0]]
    vvym = vvym[1:] + [vvym[0]]
    vvzm = vvzm[1:] + [vvzm[0]]
    vvx2 = vvx2[1:] + [vvx2[0]]
    vvy2 = vvy2[1:] + [vvy2[0]]
    vvz2 = vvz2[1:] + [vvz2[0]]
    ll = ll[1:] + [ll[0]]

    # time
    tt = []
    for i,x in enumerate(xx):
        tt.append(np.linspace(0., 1., len(x)))

    # Plot trajectories
    PlotTraj(xx, yy, ll, vf.shape)
    PlotTrajFld(tt, xx, ll, "t", "x", "trajx.pdf", ystep=1/16, vmin=0, vmax=1)
    PlotTrajFld(tt, yy, ll, "t", "y", "trajy.pdf", ystep=1/16, vmin=0, vmax=1)
    # Plot pressure jump
    PlotTrajFld(tt, ee, ll, "t", "p", "trajp.pdf", vmin=0., vmax=2.)
    # Plot error in velocity
    vmin = 1e-6
    vmax = 1e1
    PlotTrajFld(tt, vvxm, ll, "t", "vx", "trajvxm.pdf",
                vmin=vmin, vmax=vmax, ylog=True)
    PlotTrajFld(tt, vvym, ll, "t", "vy", "trajvym.pdf",
                vmin=vmin, vmax=vmax, ylog=True)
    PlotTrajFld(tt, vvx2, ll, "t", "vx", "trajvx2.pdf",
                vmin=vmin, vmax=vmax, ylog=True)
    PlotTrajFld(tt, vvy2, ll, "t", "vy", "trajvy2.pdf",
                vmin=vmin, vmax=vmax, ylog=True)

Main()
