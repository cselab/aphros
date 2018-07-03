#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re
from plotlib import *

# Figure with curvature comparison
# vf: volume fraction field
# kk: list of curvature fields
# ll: list of labels
def FigCmp(vf, kk, ll, po):
    # mesh
    x1,y1,hx,hy = GetGeom(vf.shape)
    x,y = GetMesh(x1, y1)

    # center of mass
    cx = (x * vf).sum() / vf.sum()
    cy = (y * vf).sum() / vf.sum()
    th = 1e-3
    ii = np.where((vf > th) & (vf < 1. - th))
    x = x[ii];   y = y[ii]
    dx = x - cx; dy = y - cy
    an = np.arctan2(dy, dx)
    deg = np.degrees(an)
    s = np.argsort(an)

    # exact curvature
    ane = np.linspace(-np.pi, np.pi, 200)
    dege = np.degrees(ane)
    dm,dm,dm,rx,ry = np.loadtxt('b.dat')
    ke = rx*ry / ((ry * np.cos(ane)) ** 2 + (rx * np.sin(ane)) ** 2) ** (3. / 2.)

    # average curvature
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

# Figure with volume fraction
# pt: path template
# suf: output name suffix
def FigK(pt):
    vf = ReadField2d(pt, 'u')

    k = ReadField2d(pt, 'k')
    if IsGerris(vf.shape):
        k *= -1
    # exclude points with zero curvature (for ch)
    k[np.where(k == 0.)] = np.nan

    fig, ax = PlotInitSq()
    x1,y1,hx,hy = GetGeom(vf.shape)
    xn1,yn1 = GetMeshNodes(hx, hy)
    x,y = GetMesh(x1, y1)

    # grid
    PlotGrid(ax, xn1, yn1)
    reff = GetReff(vf)
    # average curvature
    ke = 1. / reff
    print(ke)

    PlotFieldBwr(ax, k / ke, vmin=0.8, vmax=1.2)
    print(np.nanmean(k))

    # save
    po = GetFieldPath(pt, "k", "pdf")
    print(po)
    PlotSave(fig, ax, po)

# Figure with volume fraction
# pt: path template
# suf: output name suffix
def FigVf(pt):
    a = ReadField2d(pt, 'a')
    nx = ReadField2d(pt, 'nx')
    ny = ReadField2d(pt, 'ny')
    vf = ReadField2d(pt, 'u')

    fig, ax = PlotInitSq()
    x1,y1,hx,hy = GetGeom(vf.shape)
    xn1,yn1 = GetMeshNodes(hx, hy)
    x,y = GetMesh(x1, y1)

    # grid
    PlotGrid(ax, xn1, yn1)
    # field
    PlotFieldGray(ax, vf, vmin=0., vmax=1.)
    # lines
    if all([e is not None for e in [a, nx, ny]]):
        l = GetLines(x, y, a, nx, ny, hx, hy, vf)
        PlotLines(ax, *l)
    # partilces
    pa = GetFieldPath(pt, "partit")
    if os.path.isfile(pa):
        print(pa)
        PlotPart(ax, pa, sk=5)

    # save
    po = GetFieldPath(pt, "vf", "pdf")
    print(po)
    PlotSave(fig, ax, po)


def Glob(d):
    return sorted(glob.glob(os.path.join(d, "u*.dat")))

def Main():
    # directories
    dd = ['ch', 'ge']
    # labels
    ll = ['mfer', 'gerris']
    # curvature arrays
    kk = []

    for d,l in zip(dd, ll):
        print(d)
        pp = Glob(d)
        assert pp is not None
        pt = GetPathTemplate(pp[-1])

        # volume fraction
        FigVf(pt)
        # curvature
        FigK(pt)
        # title
        vf = ReadField2d(pt, "u")
        cx,cy,cz,rx,ry = np.loadtxt('b.dat')
        x1,y1,hx,hy = GetGeom(vf.shape)
        FigTitle("rx/h={:0.3f} ry/h={:0.3f}".format(rx / hx, ry / hx),
                 GetFieldPath(pt, "ttl", "pdf"))

        # append curvature
        k = ReadField2d(pt, "k")
        assert k is not None
        if IsGerris(vf.shape):
            k = k * (-1)
        kk.append(k)

    # curvature comparison
    pp = Glob(dd[0])
    pt = GetPathTemplate(pp[-1])
    vf = ReadField2d(pt, "u")
    po = 'cmp.pdf'
    FigCmp(vf, kk, ll, po)

