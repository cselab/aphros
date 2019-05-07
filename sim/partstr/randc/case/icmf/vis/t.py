#!/usr/bin/env python

import glob
import numpy as np
import os
import pickle
import math

cache = False

def GetCache(name):
    if not cache:
        return None

    c = "cache_{:}.bin".format(name)

    if os.path.isfile(c):
        print("Cache hit {:}".format(c))
        with open(c, 'rb') as f:
            return pickle.load(f)
    print("Cache miss {:}".format(c))
    return None

def SaveCache(name, r):
    if not cache:
        return r

    c = "cache_{:}.bin".format(name)

    print("Save to cache: {:}".format(c))
    with open(c, 'wb') as f:
        pickle.dump(r, f)
    return r

def Collect(f):
    cn = "partit"
    r = GetCache(cn)
    #if r is not None: return r

    u = np.genfromtxt(f, names=True, delimiter=',')
    xx = u['x']
    yy = u['y']
    zz = u['z']
    cc = u['c'].astype(int)
    kk = u['k']
    r = xx, yy, zz, cc, kk
    return SaveCache(cn, r)

# array of lines:
# line: xx, yy, zz, c
# fn: filename
def Vtk(ll, fn, comment=""):
    with open(fn, 'w') as o:
        def W(m): o.write(m + '\n')
        W("# vtk DataFile Version 2.0")
        W(comment)
        W("ASCII")
        W("DATASET POLYDATA")

        W("POINTS {:} float".format(sum([len(l[0]) for l in ll])))
        for l in ll:
            xx,yy,zz,c = l
            for i in range(len(xx)):
                W("{:} {:} {:}".format(xx[i], yy[i], zz[i]))

        W("LINES {:} {:}".format(len(ll), sum([len(l[0]) + 1 for l in ll])))
        p = 0
        for l in ll:
            W(' '.join(list(map(str, [len(l[0])] + [p + q for q in range(len(l[0]))]))))
            p += len(l[0])

        W("CELL_DATA {:}".format(len(ll)))
        W("SCALARS c float")
        W("LOOKUP_TABLE default")
        for l in ll:
            W("{:}".format(l[3]))


def IsClose(a, b):
    return abs(a - b) < 1e-5


f = sorted(glob.glob("partit*.csv"))[0]

xx, yy, zz, cc, kk = Collect(f)

ll = []
cprev = None
nnp = 9     # particles per string
for i in range(len(xx)):
    c = cc[i]
    if cprev is None or c != cprev or len(l[0]) >= nnp:
        if cprev is not None:
            ll.append(l)
        l = [[], [], [], c]
        cprev = c

    l[0].append(xx[i])
    l[1].append(yy[i])
    l[2].append(zz[i])

Vtk(ll, 'o.vtk', comment="lines connecting paritcles")
