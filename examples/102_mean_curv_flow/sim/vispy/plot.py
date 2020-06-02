#!/usr/bin/env python3

import readvtk
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def WE(msg):
    sys.stderr.write(str(msg) + "\n")

av = sys.argv

if len(av) < 2:
    sys.stderr.write('''./plot.py VTK [COLORS]
Plots connected 2D shapes in VTK as filled polygons with different colors.
VTK: path polygons from marching cubes (e.g. sm_0000.vtk)
COLORS: path to table with columns RGB, line cl defining color of index cl
Output: creates pdf (e.g. sm_0000.pdf)
''')
    exit(1)

vtk_path = av[1]
colors_path = av[2] if len(av) > 2 else ""

xx,poly,fields = readvtk.ReadVtkPoly(vtk_path)

fcl = fields['cl'].astype(int)
fclu = np.unique(fcl).astype(int)

poly = poly.astype(int)
xx = xx[:,:2]

center_cl = dict()
for cl in fclu:
    selp = np.where(fcl == cl)[0]
    if not selp.size:
        continue
    selx = poly[selp, 1:4].flatten()
    center_cl[cl] = xx[selx].mean(axis=0)

vcorn = [[0,0], [1,0], [0,1], [1,1]]
corn_cl = dict()
for corn in vcorn:
    corn = np.array(corn)
    dist_cl = {cl:sum((corn - center_cl[cl])**2) for cl in center_cl}
    cl = min(dist_cl, key=dist_cl.get)
    corn_cl[cl] = corn

resx = 256
dpi = 100
fig = plt.figure(figsize=(resx/dpi,resx/dpi))
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_axis_off()
fig.add_axes(ax)

if colors_path:
    colors = np.loadtxt(colors_path) if colors_path else []
else:
    colors = ["C" + str(int(cl)) for cl in fclu]

def GetColor(cl):
    if cl < len(colors):
        return colors[cl]
    WE("warning: fall back to C% color for cl={:}".format(cl))
    return "C" + str(int(cl))

for cl in fclu:
    selp = np.where(fcl == cl)[0]
    if not selp.size:
        continue
    selx = poly[selp, 1:4].flatten()
    xx0 = np.unique(xx[selx], axis=0)

    x = np.copy(xx0[:,0])
    y = np.copy(xx0[:,1])
    xc = np.mean(x)
    yc = np.mean(y)
    if cl in corn_cl:
        corn = corn_cl[cl]
        xc = np.mean(x)
        yc = np.mean(y)
        w = 0.5
        xc = w * xc + (1 - w) * corn[0]
        yc = w * yc + (1 - w) * corn[1]
        x = np.append(x, corn[0])
        y = np.append(y, corn[1])
    angle = np.arctan2(x - xc, y - yc)
    srt = np.argsort(angle)
    x = x[srt]
    y = y[srt]

    ax.fill(x, y, c=GetColor(cl))
    x = np.append(x, x[0])
    y = np.append(y, y[0])
    ax.plot(x, y, c='black', lw=0.8)
fig.savefig(os.path.splitext(os.path.basename(vtk_path))[0] + ".pdf")
