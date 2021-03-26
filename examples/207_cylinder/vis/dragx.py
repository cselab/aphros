#!/usr/bin/env python

try:
    import plottools
    import matplotlib.pyplot as plt
    plottools.apply_params(plt)
except:
    import matplotlib.pyplot as plt

import sys
import os
import re
import numpy as np

def getpar(name, text):
    try:
        return float(
            re.findall("set double {} ([0-9.]*)".format(name), text)[0])
    except:
        pass
    try:
        return int(re.findall("set int {} ([0-9]*)".format(name), text)[0])
    except:
        pass
    return None

def GetDrag(basedir, skip=10):
    fpath = os.path.join(basedir, "stat.dat")
    d = dict()
    ppath = os.path.join(basedir, "par.py")
    with open(ppath) as f:
        exec(f.read(), d)
    nx = d['nx']
    extent = d.get('extent', 18)
    cpd = nx / extent
    stat = np.genfromtxt(fpath, names=True)[::skip]
    return stat['t'], stat['dragx'], stat['dragy'], cpd, os.path.basename(basedir)


fig,ax = plt.subplots()

fref = "ref/dragx_ba"
t,dragx = np.genfromtxt(fref).T
ax.plot(t, dragx * 2, label="basilisk $D/h={:.0f}$".format(2**13/18), c='g')

fref = "ref/dragx_k95"
ref = np.genfromtxt(fref)
refx,refy = ref.T
refx = refx * 0.5
refy *= 1
ax.plot(refx, refy, label="koumoutsakos1995", c='0')


for path in sys.argv[1:]:
    t, dragx, dragy, cpd, name = GetDrag(path)
    ax.plot(t, dragx * 2, label="{} $D/h={:.0f}$".format(name, cpd), lw=0.75)

ax.set_xlim(0, 3)
ax.set_ylim(0, 2)
ax.set_xlabel("$t$")
ax.set_ylabel(r"$C_D$")
ax.legend()
fpath = "dragx.pdf"
print(fpath)
plt.tight_layout()
plt.savefig(fpath)
