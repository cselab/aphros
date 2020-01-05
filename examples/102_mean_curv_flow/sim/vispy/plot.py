#!/usr/bin/env python

import readvtk
import matplotlib.pyplot as plt
import numpy as np

fn = "sm_0200.vtk"

xx,poly,fields = readvtk.ReadVtkPoly(fn)

fcl = fields['cl']
fclu = np.unique(fcl)

poly = poly.astype(int)
xx = xx[:,:2]


center_cl = dict()
for cl in fclu:
    selp = np.where(fcl == cl)[0]
    if not selp.size:
        continue
    selx = poly[selp, 1]
    center_cl[cl] = xx[selx].mean(axis=0)

vcorn = [[0,0], [1,0], [0,1], [1,1]]
corn_cl = dict()
for corn in vcorn:
    corn = np.array(corn)
    dist_cl = {cl:sum((corn - center_cl[cl])**2) for cl in center_cl}
    cl = min(dist_cl, key=dist_cl.get)
    corn_cl[cl] = corn

for cl in np.unique(fcl):
    selp = np.where(fcl == cl)[0]
    if not selp.size:
        continue
    selx = poly[selp, 1:4].flatten()
    xx0 = np.unique(xx[selx], axis=0)

    x = np.copy(xx0[:,0])
    y = np.copy(xx0[:,1])
    if cl in corn_cl:
        x = np.append(x, corn_cl[cl][0])
        y = np.append(y, corn_cl[cl][1])
    xc = x.mean()
    yc = y.mean()
    angle = np.arctan2(x - xc, y - yc)
    srt = np.argsort(angle)
    x = x[srt]
    y = y[srt]

    plt.fill(x, y)
    x = np.append(x, x[0])
    y = np.append(y, y[0])
    plt.plot(x, y, c='black', lw=1)
ax = plt.gca()
ax.set_aspect('equal')
plt.savefig("a.jpg")
