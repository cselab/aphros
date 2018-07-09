#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re
from plotlib import *

# Figure with curvature depending on angle
# vf: volume fraction field
# kk: list of curvature fields
# ll: list of labels
def FigAng(vf, kk, ll, po):
    dim = GetDim(vf.shape)
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    x,y,z = GetMesh(x1, y1, z1)

    # interface cells
    th = 0
    ii = np.where((vf > th) & (vf < 1. - th))
    x = x[ii]; y = y[ii]; z = z[ii]

    # exact curvature
    ke = GetExactK(dim, x, y, z)
    # average curvature
    kea = ke.mean()

    # hires angle (high resolution)
    anh = np.linspace(-np.pi, np.pi, 200)
    degh = np.degrees(anh)
    # hires points
    cx,cy,cz,rx,ry = np.loadtxt('b.dat')
    xh = cx + np.cos(anh)
    yh = cy + np.sin(anh)
    # hires curvature
    keh = GetExactK(dim, xh, yh, cz)

    dx = x - cx ; dy = y - cy
    an = np.arctan2(dy, dx)
    deg = np.degrees(an)
    s = np.argsort(an)

    # plot exact curvature
    fig, ax = plt.subplots()
    ax.plot(degh, keh / kea, label="exact", c="0.5")

    # plot estimates
    for k,l in zip(kk,ll):
        if k is not None:
            ax.plot(deg[s], k[ii][s] / kea, label=l)

    ax.set_xlabel(r"angle [deg]")
    ax.set_ylabel(r"normalized curvature")
    q=180.
    ax.set_xticks(np.arange(-q, q*1.1, 90.))
    ax.set_xlim(-q, q)
    ax.legend()
    ax.set_ylim(0., 2.)
    ax.grid()
    plt.title("rx/h={:.3f}, ry/h={:.3f}".format(rx / hx, ry / hx))
    PlotSave(fig, ax, po)


# Figure with volume fraction
# pt: path template
# suf: output name suffix
def FigK(pt):
    vf = ReadField(pt, 'u')

    k = ReadField(pt, 'k')
    if IsGerris(vf.shape):
        k *= -1
    # exclude points with zero curvature (for ch)
    k[np.where(k == 0.)] = np.nan

    fig, ax = PlotInitSq()
    # q: dummy
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    xn1,yn1,zn1 = GetMeshNodes(hx, hy, hz)
    x,y,z = GetMesh(x1, y1, z1)

    # grid
    PlotGrid(ax, xn1, yn1)
    reff = GetReff(vf, GetDim(k.shape))
    # average curvature
    ke = 1. / reff

    PlotFieldBwr(ax, k / ke, vmin=0.8, vmax=1.2)

    # save
    po = GetFieldPath(pt, "k", "pdf")
    print(po)
    PlotSave(fig, ax, po)

# Figure with volume fraction
# pt: path template
# suf: output name suffix
def FigVf(pt):
    a = ReadField(pt, 'a')
    nx = ReadField(pt, 'nx')
    ny = ReadField(pt, 'ny')
    vf = ReadField(pt, 'u')

    fig, ax = PlotInitSq()
    # q: dummy
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    xn1,yn1,zn1 = GetMeshNodes(hx, hy, hz)
    x,y,z = GetMesh(x1, y1, z1)

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
        PlotPart(ax, pa, sk=4)

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

        # append curvature
        k = ReadField(pt, "k")
        assert k is not None
        if IsGerris(k.shape):
            k = k * (-1)
        kk.append(k)

    # volume fraction from dd[0]
    pp = Glob(dd[0])
    pt = GetPathTemplate(pp[-1])
    vf = ReadField(pt, "u")

    # title
    cx,cy,cz,rx,ry = np.loadtxt('b.dat')
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    po = 'ttl.pdf'
    FigTitle("rx/h={:0.3f} ry/h={:0.3f}".format(rx / hx, ry / hx), po)

    # curvature vs angle
    po = 'ang.pdf'
    FigAng(vf, kk, ll, po)

    # curvature histogram
    po = 'hist.pdf'
    FigHistK(vf, kk, ll, po)

Main()
