#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import pickle
import glob
import os
import scipy.optimize

def GetCache(name):
    c = "cache_{:}.bin".format(name)
    if os.path.isfile(c):
        print("Cache hit {:}".format(c))
        with open(c, 'rb') as f:
            return pickle.load(f)
    print("Cache miss {:}".format(c))
    return None

def SaveCache(name, r):
    c = "cache_{:}.bin".format(name)
    print("Save to cache: {:}".format(c))
    with open(c, 'wb') as f:
        pickle.dump(r, f)
    return r

def Read(p):
    return np.genfromtxt(p, delimiter=',', names=True)

def GetMap(d):
    if d.size == 1:
        d = np.array([d])
    return {row['cl']:row for row in d}

def DiffRel(r0, r1):
    return (r0 - r1) / max(r0, r1) > 0.1

# True if reduction of volume is less 50%
def VolHalf(v0, v1):
    return (v0 - v1) / max(v0, v1) < 0.5

def Write(d, ccl, i, p="a_%.csv"):
    m = GetMap(d)
    dc =  np.array([m[cl] for cl in ccl], dtype=d.dtype)
    p = p.replace('%', "{:04d}".format(i))
    print(p)
    np.savetxt(p, dc, header=','.join(d.dtype.names), comments='', delimiter=',')

def StrHeader(d):
    return ','.join(d.dtype.names)

def StrRow(row):
    return ','.join(map(str, row))

def WriteLine(s, out):
    out.write(s + '\n')

def Get(m, ccl, field):
    return np.array([m[cl][field] for cl in ccl if cl in m])

def ApplyRemap(ccl, rmp):
    return [rmp.get(ccl[i], ccl[i]) for cl in ccl]

def GetGood(p0, p1):
    print(p0, p1)
    d0 = Read(p0)
    d1 = Read(p1)

    m0 = GetMap(d0)
    m1 = GetMap(d1)

    ccl0 = set(m0.keys())
    ccl1 = set(m1.keys())

    dis = ccl0 - ccl1
    new = ccl1 - ccl0
    diff = {cl for cl in ccl0 & ccl1 if DiffRel(m0[cl]['vf'],  m1[cl]['vf'])}
    good = ccl0 & ccl1 - diff
    bad = diff | new
    Plot(m0, m1, bad)
    Plot(m0, m1, dis)
    rmp = dict()
    def Cost(c0, c1):
        pass
    print(dis)
    print(new)

    #Plot(m0, m1, dis.union(new))
    exit()
    return good, rmp

def Plot(m0, m1, diff):
    v='x'
    vv='y'
    plt.scatter(
            Get(m1, diff, v),
            Get(m1, diff, vv),
            #Get(m1, diff, 'cl'),
            #Get(m1, diff, 'y'),
            s=Get(m1, diff, 'r')**2*2e5)
    plt.scatter(
            #Get(m0, diff, v) + Get(m0, diff, 'v'+v)*dt,
            #Get(m0, diff, vv) + Get(m0, diff, 'v'+vv)*dt,
            Get(m0, diff, v),
            Get(m0, diff, vv),
            #Get(m0, diff, 'cl'),
            #Get(m0, diff, 'y') + Get(m0, diff, 'vy')*dt,
            s=Get(m0, diff, 'r')**2*2e5, edgecolor='red', facecolor='none')
    #plt.scatter(Get(m0, dis, 'x'), Get(m0, dis, 'y'))
    #plt.scatter(Get(m1, new, 'z'), Get(m1, new, 'y'), s=15)
    plt.savefig("a.pdf")

# p: path traj_0000.csv
def GetIndex(path):
    return int(re.findall(".*_(\d*)", path)[0])

cl0 = float(sys.argv[1])
pp = sys.argv[2:]

out = sys.stdout

cl = float(clinit)
m1 = None
m0 = None

WriteLine('n,' + StrHeader(Read(pp[0])), out)

for i in range(0, len(pp)):
    print(i)
    m0 = m1
    m1 = GetMap(Read(pp[i]))
    if cl in m1 and (m0 is None or cl not in m0
            or VolHalf(m0[cl]['vf'], m1[cl]['vf'])):
        WriteLine(str(i) + ',' + StrRow(m1[cl]), out)
    else:
        if m0 is None or cl not in m0:
            continue
        def Cost(cln):
            dist = ((m0[cl]['x'] - m1[cln]['x'])**2 +
                  (m0[cl]['y'] - m1[cln]['y'])**2 +
                  (m0[cl]['z'] - m1[cln]['z'])**2) ** 0.5
            r0 = m0[cl]['r']
            r1 = m1[cln]['r']
            return dist + min(0, r0 - r1)
        new = set(m1) - set(m0)
        if not len(new):
            continue
        clbest = list(new)[0]
        for cln in new:
            if Cost(cln) < Cost(clbest):
                clbest = cln
        cl = clbest
        print('cl=', cl)
        WriteLine(str(i) + ',' + StrRow(m1[cl]), out)

