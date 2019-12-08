#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import pickle
import glob
import os
import scipy.optimize

dt = 0.05
nx = 96
h = 1 / nx

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
    return {d['cl'][i]:d[i] for i in range(0, len(d)) if d['r'][i] > h}

def DiffRel(r0, r1):
    return (r0 - r1) / max(r0, r1) > 0.1

def Write(d, ccl, i, p="a_%.csv"):
    m = GetMap(d)
    dc =  np.array([m[cl] for cl in ccl], dtype=d.dtype)
    p = p.replace('%', "{:04d}".format(i))
    print(p)
    np.savetxt(p, dc, header=','.join(d.dtype.names), comments='', delimiter=',')

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

#av = sys.argv[1:]
av = sorted(glob.glob("traj_00??.csv"))

namecache = "goodall"
goodall = GetCache(namecache)
goodall = None

if goodall is None:
    goodall = set(GetMap(Read(av[-1])).keys())
    for i in reversed(range(0, len(av) - 1)):
        p0 = av[i]
        p1 = av[i + 1]
        good = GetGood(p0, p1)
        Write(Read(p1), goodall - good, i + 1, "b_%.csv")
        goodall = goodall.intersection(good)
        Write(Read(p0), goodall, i)
    SaveCache(namecache, goodall)


'''
for i in range(0, len(av)):
    p = av[i]
    print(p)
    i = GetIndex(p)
    Write(Read(p), goodall, i)
'''
