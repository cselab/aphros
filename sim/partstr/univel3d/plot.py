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

# Figure with volume fraction
# pt: path template
# suf: output name suffix
def FigVf(pt):
    vf = ReadField(pt, 'vf')

    # center of mass
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    x,y,z = GetMesh(x1, y1, z1)
    ii = np.where(np.isfinite(vf))
    cx = np.sum(x[ii] * vf[ii]) / np.sum(vf[ii])
    cy = np.sum(y[ii] * vf[ii]) / np.sum(vf[ii])
    cz = np.sum(z[ii] * vf[ii]) / np.sum(vf[ii])
    iz = np.argmin(abs(cz - z))
    print("FigVf: c={:}, ciz={:}".format((cx, cy, cz), iz))

    # slice through center of mass
    if vf is not None: vf = vf[:,:,iz:iz+1]

    fig, ax = PlotInitSq()
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    xn1,yn1,zn1 = GetMeshNodes(hx, hy, hz)
    x,y,z = GetMesh(x1, y1, z1)

    # grid
    PlotGrid(ax, xn1, yn1)
    # field
    PlotFieldGray(ax, vf, vmin=0., vmax=1.)
    # save
    po = GetFieldPath(pt, "vf", "pdf")
    print(po)
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
    cx = [] ; cy = [] ; cz = []
    for i,p in enumerate(pp):
        # report
        sys.stdout.write("{:} ".format(i))
        sys.stdout.flush()
        # read array
        vf = ReadArray(p)
        x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
        x,y,z = GetMesh(x1, y1, z1)
        # compute center
        cx0 = (x * vf).sum() / vf.sum()
        cy0 = (y * vf).sum() / vf.sum()
        cz0 = (z * vf).sum() / vf.sum()
        cx.append(cx0)
        cy.append(cy0)
        cz.append(cz0)
    sys.stdout.write("\n")
    sys.stdout.flush()

    cx = np.array(cx)
    cy = np.array(cy)
    cz = np.array(cz)
    return cx,cy,cz

# Extracts value of field along trajectory.
# pp: list of paths to volume fraction fields
# cx,cy: trajectory
# fld: field prefix (e.g. 'p')
# Returns:
# cu: field at center
def GetTrajFld(pp, cx, cy, cz, fld):
    cu = []
    for i,p in enumerate(pp):
        # report
        sys.stdout.write("{:} ".format(i))
        sys.stdout.flush()
        # path template
        pt = GetPathTemplate(p)
        # read field
        u = ReadField(pt, fld)
        # mesh
        x1,y1,z1,hx,hy,hz = GetGeom(u.shape)
        x,y,z = GetMesh(x1, y1, z1)
        # index
        ci0 = int(cx[i] / hx)
        cj0 = int(cy[i] / hy)
        ck0 = int(cz[i] / hy)
        # field
        cu.append(u[ci0, cj0, ck0])
    sys.stdout.write("\n")
    sys.stdout.flush()

    cu = np.array(cu)
    return cu

# Plots trajectories
# xx,yy: list of arrays for axes
# ll: line labels
# lx,ly: axes labels
# po: output path
def PlotTrajFld(xx, yy, ll, lx, ly, po, vmin=None, vmax=None, ystep=None):
    global reff, hx
    fig, ax = PlotInit()
    if vmin is None: vmin = yy[0].min()
    if vmax is None: vmax = yy[0].max()
    i = 0
    for x,y,l in zip(xx, yy, ll):
        if i == len(xx) - 1: # separate for last
            ax.plot(x, y, label=l, c="0.5", ls='--')
        else:
            ax.plot(x, y, label=l)
        i += 1
    ax.legend()
    ax.grid(True)
    ax.set_xlabel(lx)
    ax.set_ylabel(ly)
    ax.set_ylim(vmin, vmax)
    if ystep is not None:
        ax.set_yticks(np.arange(vmin, vmax + 1e-10, ystep))
    PlotSave(fig, ax, po)

# Plots trajectories
# xx,yy: list of arrays
# ll: labels
# s: shape of of data array
def PlotTraj(xx, yy, ll, s):
    global reff, hx
    fig, ax = PlotInitSq()
    ax.set_aspect('equal')
    i = 0
    for x,y,l in zip(xx, yy, ll):
        if i == len(xx) - 1: # separate for last
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
    # pressure along trajectories
    ee = []

    # volume fraction from dd[0]
    pp = Glob(dd[0], "vf")
    pt = GetPathTemplate(pp[0])
    vf = ReadField(pt, "vf")

    cx,cy,cz,rx,ry,rz = LoadBub()
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)

    # exact trajectories
    x,y = GetTrajE(0.3, 0.3, 0.4, 0.3, 1.)
    # exact pressure jump
    sig = 1.
    eex = 2. * sig / rx
    e = x * 0 + eex
    # relative to exact
    e /= eex
    # append
    xx.append(x)
    yy.append(y)
    ee.append(e)
    ll.append("exact")

    # lines
    ll.append('mfer')
    ll.append('gerris')
    for d,l in zip(dd, ll):
        print(d)
        pp = Glob(d, "vf")
        assert len(pp)
        pt = GetPathTemplate(pp[0])
        # volume fraction
        FigVf(pt)
        # trajectory
        x,y,z = GetTraj(pp)
        # pressure
        e = GetTrajFld(pp, x, y, z, "p")
        # pressure at corner
        e0 = GetTrajFld(pp, x * 0, y * 0, z * 0, "p")
        # pressure jump
        e -= e0
        # relative to exact
        e /= eex
        xx.append(x[1:])
        yy.append(y[1:])
        ee.append(e[1:])

    # reorder lines
    xx = xx[1:] + [xx[0]]
    yy = yy[1:] + [yy[0]]
    ee = ee[1:] + [ee[0]]
    ll = ll[1:] + [ll[0]]

    # time
    tt = []
    for i,x in enumerate(xx):
        tt.append(np.linspace(0., 1., len(x)))


    # Plot trajectories
    PlotTraj(xx, yy, ll, vf.shape)
    # Plot pressure jump
    PlotTrajFld(tt, ee, ll, "t", "p", "trajp.pdf", vmin=0., vmax=2.)
    # Plot x,y
    PlotTrajFld(tt, xx, ll, "t", "x", "trajx.pdf", ystep=0.1, vmin=0.25, vmax=0.75)
    PlotTrajFld(tt, yy, ll, "t", "y", "trajy.pdf", ystep=0.1, vmin=0.25, vmax=0.65)

Main()
