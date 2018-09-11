#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob
import os
import re
from ch.plot import *
import scipy
import scipy.interpolate


# Figure with curvature depending on angle
# vf: volume fraction field
# kk: list of curvature fields of shape vf.shape
# ll: list of labels
def FigAng(vf, kk, ll, po):
    dim = GetDim(vf.shape)
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    x,y,z = GetMesh(x1, y1, z1)

    # interface cells
    th = 0
    ii = np.where((vf > th) & (vf < 1. - th))
    x = x[ii]; y = y[ii]; z = z[ii]

    # hires angle (high resolution)
    uh = np.linspace(-np.pi, np.pi, 200)
    vh = np.linspace(-np.pi * 0.5, np.pi * 0.5, 200)
    #uh, vh = np.meshgrid(uh, vh, indexing='ij')
    deguh = np.degrees(uh)
    degvh = np.degrees(vh)
    # hires points
    # spherical coordinates
    cx,cy,cz,rx,ry,rz = LoadBub()
    r = ((np.cos(uh) * np.cos(vh) / rx) ** 2 +
         (np.sin(uh) * np.cos(vh) / ry) ** 2 +
         (np.sin(vh) / rz) ** 2) ** (-0.5)
    xh = cx + r * np.cos(uh) * np.cos(vh)
    yh = cy + r * np.sin(uh) * np.cos(vh)
    zh = cz + r * np.sin(vh)
    # hires curvature
    keh = GetExactK(dim, xh.flatten(), yh.flatten(), zh.flatten())
    keh = keh.reshape(xh.shape)
    # average curvature
    kea = keh.mean()

    # cartesian to spherical
    dx = x - cx ; dy = y - cy ; dz = z - cz
    u = np.arctan2(dy, dx)
    degu = np.degrees(u)
    v = np.arctan2(dz, (dx ** 2 + dy ** 2) ** 0.5)
    degv = np.degrees(v)

    # line slice
    us = np.linspace(-np.pi, np.pi, 200)
    vs = np.linspace(-np.pi * 0.5, np.pi * 0.5, 200)
    # deg
    degus = np.degrees(us)
    degvs = np.degrees(vs)
    kks = []
    for k in kk:
        if k is not None:
            kks.append(scipy.interpolate.griddata(
                    (degu, degv), k[ii], (degus, degvs)))
        else:
            kks.append(None)

    fig, ax = PlotInit()

    # plot estimates
    for k,l in zip(kks,ll):
        if k is not None:
            ax.plot(degus, k / kea, label=l)

    # plot lowres exact curvature
    ke = GetExactK(dim, x, y, z)
    kes = scipy.interpolate.griddata((degu, degv), ke, (degus, degvs))
    ax.plot(degus, kes / kea, label="exact")

    # plot exact curvature
    ax.plot(degus, keh / kea, label="exact", c="0.5", ls='--')

    ax.set_xlabel(r"$\varphi$ [deg]")
    ax.set_ylabel(r"normalized curvature")
    ax.set_ylim(0., 2.)
    ax.grid()
    ax.legend()
    #plt.show()
    PlotSave(fig, ax, po)


# Figure with volume fraction
# pt: path template
# suf: output name suffix
def FigVf(pt):
    vf = ReadField(pt, 'u')

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
    cx,cy,cz,rx,ry,rz = LoadBub()
    x1,y1,z1,hx,hy,hz = GetGeom(vf.shape)
    title = "rx/h={:0.3f} ry/h={:0.3f} rz/h={:0.3f}".format(
            rx / hx, ry / hx, rz / hx)

    # curvature vs angle
    po = 'ang.pdf'
    FigAng(vf, kk, ll, po)

    # curvature histogram
    po = 'hist.pdf'
    FigHistK(vf, kk, ll, po, title=title)

Main()
