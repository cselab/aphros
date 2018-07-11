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
        print(pa)
        PlotPart(ax, pa, sk=4)

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
    cx = [] ; cy = []
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
        cx.append(cx0)
        cy.append(cy0)
    sys.stdout.write("\n")
    sys.stdout.flush()

    cx = np.array(cx)
    cy = np.array(cy)
    return cx,cy

# Extracts value of field along trajectory.
# pp: list of paths to volume fraction fields
# cx,cy: trajectory
# fld: field prefix (e.g. 'p')
# Returns:
# cu: field at center
def GetTrajFld(pp, cx, cy, fld):
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
        # field
        cu.append(u[ci0, cj0, 0])
    sys.stdout.write("\n")
    sys.stdout.flush()

    cu = np.array(cu)
    return cu

# Plots trajectories
# xx,yy: list of arrays for axes
# ll: line labels
# lx,ly: axes labels
# po: output path
def PlotTrajFld(xx, yy, ll, lx, ly, po):
    global reff, hx
    fig, ax = plt.subplots()
    for x,y,l in zip(xx, yy, ll):
        ax.plot(x, y, label=l)
    ax.legend()
    ax.set_xlabel(lx)
    ax.set_ylabel(ly)
    PlotSave(fig, ax, po)

# Plots trajectories
# xx,yy: list of arrays
# ll: labels
# s: shape of of data array
def PlotTraj(xx, yy, ll, s):
    global reff, hx
    fig, ax = PlotInitSq()
    ax.set_aspect('equal')
    for x,y,l in zip(xx, yy, ll):
        ax.plot(x, y, label=l)
    ax.legend()
    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    x1,y1,z1,hx,hy,hz = GetGeom(s)
    xn1,yn1,zn1 = GetMeshNodes(hx, hy, hz)
    PlotGrid(ax, xn1, yn1)

    q,q,q,rx,ry = np.loadtxt('b.dat')
    x1,y1,z1,hx,hy,hz = GetGeom(s)
    po = 'ttl.pdf'
    FigTitle("r/h={:0.3f}".format(rx / hx), po)
    plt.title("r/h={:0.3f}".format(rx / hx))
    po = 'traj.pdf'
    PlotSave(fig, ax, po)

def Main():
    # directories
    dd = ['ch', 'ge']
    # labels
    ll = ['mfer', 'gerris']
    # trajectories
    xx = []
    yy = []
    # pressure along trajectories
    ee = []

    for d,l in zip(dd, ll):
        print(d)
        pp = Glob(d, "vf")
        assert len(pp)
        pt = GetPathTemplate(pp[0])
        # volume fraction
        FigVf(pt)
        # trajectory
        x,y = GetTraj(pp)
        # pressure along
        e = GetTrajFld(pp, x, y, "p")
        # pressure at corner
        e0 = GetTrajFld(pp, x * 0, y * 0, "p")
        # pressure jump along
        e -= e0
        xx.append(x)
        yy.append(y)
        ee.append(e)

    # volume fraction from dd[0]
    pp = Glob(dd[0], "vf")
    pt = GetPathTemplate(pp[0])
    vf = ReadField(pt, "vf")

    # title
    cx,cy,cz,rx,ry = np.loadtxt('b.dat')
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    po = 'ttl.pdf'
    FigTitle("r/h={:0.3f}".format(rx / hx), po)

    # exact trajectories
    x,y = GetTrajE(0.3, 0.3, 0.4, 0.3, 1.)
    # exact pressure jump
    sig = 1.
    e = x * 0 + sig / rx
    # append
    xx = [x] + xx
    yy = [y] + yy
    ee = [e] + ee
    ll = ["exact"] + ll

    # time
    tt = []
    for i,x in enumerate(xx):
        tt.append(np.linspace(0., 1., len(x)))

    # Plot trajectories
    PlotTraj(xx, yy, ll, vf.shape)
    # Plot pressure along
    PlotTrajFld(tt, ee, ll, "t", "p", "trajp.pdf")
    # Plot x,y along
    PlotTrajFld(tt, xx, ll, "t", "x", "trajx.pdf")
    PlotTrajFld(tt, yy, ll, "t", "y", "trajy.pdf")

Main()
