#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re
from plotlib import *


# Figure with volume fraction
# pt: path template
# suf: output name suffix
def FigVf(pt):
    # load 3d field
    vf = ReadField3d(pt, 'u')

    # center of mass
    x1,y1,z1,hx,hy,hz = GetGeom3(vf.shape)
    x,y,z = GetMesh3(x1, y1, z1)
    ii = np.where(np.isfinite(vf))
    cx = np.sum(x[ii] * vf[ii]) / np.sum(vf[ii])
    cy = np.sum(y[ii] * vf[ii]) / np.sum(vf[ii])
    cz = np.sum(z[ii] * vf[ii]) / np.sum(vf[ii])
    iz = np.argmin(abs(cz - z))
    print("FigVf: c={:}, ciz={:}".format((cx, cy, cz), iz))

    # slice through center of mass
    if vf is not None: vf = vf[:,:,iz]

    fig, ax = PlotInitSq()
    x1,y1,hx,hy = GetGeom(vf.shape)
    xn1,yn1 = GetMeshNodes(hx, hy)
    x,y = GetMesh(x1, y1)

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
        k = ReadField3d(pt, "k")
        assert k is not None
        if IsGerris(k.shape):
            k = k * (-1)
        kk.append(k)

    # volume fraction from dd[0]
    pp = Glob(dd[0])
    pt = GetPathTemplate(pp[-1])
    vf = ReadField3d(pt, "u")

    # title
    cx,cy,cz,rx,ry = np.loadtxt('b.dat')
    x1,y1,z1,hx,hy,hz = GetGeom3(vf.shape)
    po = 'ttl.pdf'
    FigTitle("rx/h={:0.3f} ry/h={:0.3f}".format(rx / hx, ry / hx), po)

    # curvature depending on angle
    po = 'hist.pdf'
    FigHistK(vf, kk, ll, po)

Main()
