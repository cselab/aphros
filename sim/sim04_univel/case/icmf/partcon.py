#!/usr/bin/env python

import glob
import numpy as np
import os
import pickle
import math
import sys

def Collect(f):
    u = np.genfromtxt(f, names=True, delimiter=',')
    xx = u['x']
    yy = u['y']
    zz = u['z']
    cc = u['c'].astype(int)
    kk = u['k']
    return xx, yy, zz, cc, kk

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

av = sys.argv
if len(av) < 3:
    sys.stderr.write('''usage: {:} src out
src: file of partit*.csv
out: output vtk
'''.format(av[0]))
    exit(1)

f = av[1]
out = av[2]

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

Vtk(ll, out, comment="lines connecting paritcles")
