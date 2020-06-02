#!/usr/bin/env python3

import numpy as np
import cv2
import os
import matplotlib.pyplot as plt
import sys
import colorsys
import scipy.stats

def kmeans(img):
    nx,ny,nc = img.shape
    Z = img.reshape((-1,3))
    Z = np.float32(Z)
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 1.0)
    K = 100  # from paper
    comp,label,center=cv2.kmeans(Z,K,None,criteria,1,cv2.KMEANS_RANDOM_CENTERS)
    center0 = center
    center = np.uint8(center)
    res = center[label.flatten()]
    res = res.reshape((img.shape))
    resl = label.flatten()
    resl = resl.reshape((ny,nx))
    return resl

def recolor(u):
    # initial color
    cl = np.arange(nx * ny).reshape(nx, ny)
    # u: image with 3 channels, shape (nx,ny,3)
    # cl: color, shape (nx,ny)
    def Iter(u, cl):
        for d in [0, 1]:
            for s in [-1, 1]:
                u_d = np.roll(u, s, axis=d)
                cl_d = np.roll(cl, s, axis=d)
                sel = np.where(np.mean(abs(u - u_d), axis=2) < 50)
                changes = np.where(cl[sel] != cl_d[sel])
                cl[sel] = np.minimum(cl[sel], cl_d[sel])
        return cl, changes[0].size

    for i in range(100):
        cl, ch = Iter(u, cl)
        print(i, ch)

    clu = np.unique(cl)
    num,bins = np.histogram(cl, bins=clu)
    sel = np.where(num > 200)
    plt.plot(num[sel])
    print(bins[sel], num[sel])
    o = "a.pdf"
    plt.savefig(o)
    plt.close()

def pool(u):
    st = 4
    nx,ny = u.shape
    uu = np.zeros((nx, ny, (st*2+1)**2))
    i = 0;
    for dx in range(-st, st+1):
        for dy in range(-st, st+1):
            uu[:,:,i] = np.pad(u, ((st, st), (st, st)), mode='reflect')[
                    st-dx:nx+st-dx, st-dy:ny+st-dy]
            i += 1
    u = scipy.stats.mode(uu, axis=2)[0]
    u = u.reshape(nx, ny)
    return u

# Write uniform grid data
# u -- 2d or 3d array
# Format:
# <Nx> <Ny> <Nz>
# <u[0,0,0]> <u[1,0,0]> ...
def WritePlain(u, path):
    s = u.shape
    assert len(s) in [2, 3]
    if (len(s) == 2):
        u = u.reshape((s[0], s[1], 1))
    with open(path, 'w') as f:
        f.write("{:} {:} {:}\n".format(*u.shape))
        u = u.flatten()
        np.savetxt(f, u, newline='', fmt='%.16g ')


fn = "cl.npy"
if True or not os.path.exists(fn):
    f = "a.png"
    u = cv2.imread(f)
    u = cv2.blur(u,(2,2))
    cl = kmeans(u)
    np.save(fn, cl)

cl = np.load(fn)
cl = pool(cl)
WritePlain(cl, "cl.dat")
plt.imshow(cl)
o = "cl.pdf"
plt.savefig(o)
plt.close()
