#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
import re

import aphros


# u: numpy array (2d or 3d slice)
def Get2d(u):
    if u is None:
        return None
    s = u.shape
    if len(s) == 2:
        return u
    else:
        assert len(s) >= 3
        return u[:, :, 0].reshape((s[0], s[1]))

def Read2d(p):
    return Get2d(aphros.plot.ReadArray(p))

pre = 'u'
pp = sorted(glob.glob(pre + "*.dat"))
#pp = pp[1:2]
for p in pp:
    suf = re.findall(pre + "(.*)", p)[0]
    print("p={:}, suf={:}".format(p, suf))

    u = Read2d(p)

    [sx, sy] = u.shape  # cells

    # check if gerris (odd sx)
    ge = (sx % 2 == 1)

    if ge:
        sx -= 1
        sy -= 1

    hx = 1. / sx
    hy = 1. / sy

    # equivalent radius
    reff = ((u.sum() * hx * hy) / np.pi)**0.5
    print("reff={:}".format(reff))

    # nodes, 1 means 1d
    xn1 = np.arange(sx + 1) * hx
    yn1 = np.arange(sy + 1) * hy
    # cell centers
    x1 = (0.5 + np.arange(sx)) * hx
    y1 = (0.5 + np.arange(sy)) * hy
    x, y = np.meshgrid(x1, y1)
    x = x.T
    y = y.T

    a = Read2d('a' + suf)  # alpha
    nx = Read2d('nx' + suf)  # normal
    ny = Read2d('ny' + suf)

    # volume fraction u
    po = os.path.splitext(p)[0] + ".pdf"
    print(po)
    fig, ax = aphros.plot.PlotInitSq()
    aphros.plot.PlotGrid(ax, xn1, yn1)
    vmax = abs(u).max()
    aphros.plot.PlotFieldGray(ax, u, vmin=0., vmax=1.)
    if all([e is not None for e in [a, nx, ny]]):
        l = aphros.plot.GetLines(x, y, a, nx, ny, hx, hy, u)
        aphros.plot.PlotLines(ax, *l)
    # plot partilces
    n = re.findall("_(\d*)\.", suf)
    if n:
        pa = "partit.{:}.csv".format(n[0])
        if os.path.isfile(pa):
            print(pa)
            PlotPart(ax, pa, sk=5)
    # save
    aphros.plot.PlotSave(fig, ax, po)
