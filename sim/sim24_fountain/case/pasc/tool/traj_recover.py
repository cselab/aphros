#!/usr/bin/env python

import numpy as np
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

clinit = float(sys.argv[1])
pp = sys.argv[2:]

out = sys.stdout

cl = float(clinit)
m1 = None
m0 = None

WriteLine('n,' + StrHeader(Read(pp[0])), out)

for i in range(0, len(pp)):
    sys.stderr.write("{:}".format(i) + '\n')
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
        WriteLine(str(i) + ',' + StrRow(m1[cl]), out)

