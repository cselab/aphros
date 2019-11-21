#!/usr/bin/env python

import glob
import numpy as np
import os
import pickle
import math

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

def GetNb(f):
    return np.genfromtxt(f, names=True, delimiter=',').shape[0]

def Collect(ff):
    cn = "traj"
    r = GetCache(cn)
    if r is not None: return r

    nb = GetNb(ff[0])
    nt = len(ff)

    xx = np.zeros((nt, nb))
    yy = np.zeros((nt, nb))
    zz = np.zeros((nt, nb))
    aa = np.zeros((nt, nb))
    rr = np.zeros((nt, nb))

    for t,f in enumerate(ff):
        u = np.genfromtxt(f, names=True, delimiter=',')
        for q in u:
            c = int(q['cl'])
            xx[t,c] = q['x']
            yy[t,c] = q['y']
            zz[t,c] = q['z']
            aa[t,c] = q['vf']
            rr[t,c] = q['r']
    r = xx, yy, zz, aa, rr
    return SaveCache(cn, r)


ff = sorted(glob.glob("traj/*.csv"))

ff = ff[::5]

xx, yy, zz, aa, rr = Collect(ff)

def GetDist(t, b, bo):
    global xx, yy, zz
    return ((xx[t, b] - xx[t, bo]) ** 2 +
            (yy[t, b] - yy[t, bo]) ** 2 +
            (zz[t, b] - zz[t, bo]) ** 2) ** 0.5

def GetCoalPartners():
    cn = "coalpart"
    r = GetCache(cn)
    if r is not None: return r

    global nb, nt, ee

    # coalescence partners, shape (nb)
    pp = np.zeros((nb)).astype(int) * 0 - 1

    # fill pp
    for b in range(nb):
        if (ee[b] < nt - 1):
            t = ee[b]
            mbo = 0
            for bo in range(nb):
                if GetDist(t, b, bo) < GetDist(t, b, mbo) and ee[bo] > ee[b] and b != bo:
                    mbo = bo
            pp[b] = mbo
        assert pp[b] != b
    r = pp
    return SaveCache(cn, r)

nt,nb = xx.shape

# lifetime, max timestep while exists , shape (nb)
ee = np.zeros((nb)).astype(int)

# vf threshold
def Ath():
    return 0. * aa[0,0]


for b in range(nb):
    ee[b] = np.max(np.where(aa[:,b] > Ath()))

# number of bubbles over time, shape (nt)
nn = np.zeros((nt))
tt = np.arange(nt) * 20 / (nt - 1)

for t in range(nt):
    nn[t] = len(np.where(ee >= t)[0])

np.savetxt("nb.dat", np.array((tt, nn)).T, header="t nb", comments='', delimiter=' ')

pp = GetCoalPartners()

# bubble indices
#sel0 = range(0, nb, 50)
sel0 = range(0,nb)
cx = np.pi * 0.5
cy = np.pi * 0.5
cz = np.pi
sel0 = np.argsort((xx[0,:]-cx)**2 +(yy[0,:]-cy)**2 +  (zz[0,:]-cz)**2)[:]
sel0 = [91]
print("selection: {:}".format(sel0))


sel = set()
# add all references as coal parnters
def S(b):
    global sel
    if b not in sel:
        sel.add(b)
        if ee[b] < nt - 1:  # if merges into other
            assert pp[b] != -1
            S(pp[b])
        for bo in range(nb): # if other merged into
            if pp[bo] == b:
                S(bo)

for b in sel0:
    S(b)

#sel = list(sel)[:10]

print("selection with coalparners: {:}".format(sel))


# write header:
def WH(o):
    o.write("cl,t,x,y,z,vf,r,ncoal\n")

def GetNcoal(a):
    return max(int(a / aa[0,0] + 0.5) - 1, 0)

# write entry
def W(t, b, o):
    o.write("{:},{:},{:},{:},{:},{:},{:},{:}\n".format(
        b, t, xx[t,b], yy[t,b], zz[t,b], aa[t,b],rr[t,b],GetNcoal(aa[t,b])))

L = 2. * math.pi

def Clip(x):
    if x < 0.:
        x += L
        Clip(x)
    elif x > L:
        x -= L
        Clip(x)
    return x


# write entry
def WW(t, b, o):
    x = Clip(xx[t,b])
    y = Clip(yy[t,b])
    z = Clip(zz[t,b])

    o.write("{:},{:},{:},{:},{:},{:},{:},{:}\n".format(
        b, t, x, y, z, aa[t,b],rr[t,b],GetNcoal(aa[t,b])))

# array of lines:
# line: xx, yy, zz, c
# fn: filename
def Vtk(ll, fn):
    with open(fn, 'w') as o:
        def W(m): o.write(m + '\n')
        W("# vtk DataFile Version 2.0")
        W("bubble trajectories")
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
        W("SCALARS ncoal float")
        W("LOOKUP_TABLE default")
        for l in ll:
            W("{:}".format(l[3]))


# sel init
with open("o0.csv", 'w') as o:
    WH(o)
    for b in sel:
        W(0, b, o)


# sel coal
with open("ocoal.csv", 'w') as o:
    WH(o)
    for b in sel:
        W(ee[b], b, o)
        bo = pp[b]
        if bo != -1:
            W(ee[b], bo, o)
            #dt = len(ff) // 30
            dt = 1
            if ee[b] + dt <= ee[bo]:
                W(ee[b] + dt, bo, o)

        if ee[b] == nt - 1:
            W(ee[b], b, o)

# sel trajectory
with open("o.csv", 'w') as o:
    WH(o)
    for b in sel:
        for t in range(nt):
            if aa[t,b] > Ath():
                W(t, b, o)


ll = []
for b in sel:
    t = 0
    while t < ee[b]:
        a = aa[t,b]
        l = [[],[],[],GetNcoal(a)]
        while t < ee[b] and GetNcoal(a) == GetNcoal(aa[t,b]):
            l[0].append(xx[t,b])
            l[1].append(yy[t,b])
            l[2].append(zz[t,b])
            t += 1
        ll.append(l)

Vtk(ll, 'o.vtk')

z0 = np.pi * 0.5
z1 = np.pi * 1.5
L = 2. * np.pi

# sel trajectory
with open("all0.csv", 'w') as o:
    WH(o)
    t = 0
    for b in range(nb):
        z = zz[t,b]
        z = Clip(z)
        if z >= z0 and z <= z1:
            WW(t, b, o)

# sel trajectory
with open("all1.csv", 'w') as o:
    WH(o)
    t = nt - 1
    for b in range(nb):
        z = zz[t,b]
        z = Clip(z)
        if z >= z0 and z <= z1:
            WW(t, b, o)

exit()

import matplotlib.pyplot as plt

plt.plot(nn)
plt.show()
