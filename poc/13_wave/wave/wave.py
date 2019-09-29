#!/usr/bin/env python

import numpy as np
import scipy
import scipy.interpolate
import scipy.sparse as sp
import scipy.sparse.linalg
import matplotlib.pyplot as plt

# solve Posson equation for
# laplace p = -vort(u,v)
# stream function p (psi) defined via u = dp/dy, v = -dp/dx
# ADHOC: periodic in x
def stream_direct(u, v):
    print("stream_direct")
    f = -vort(u,v)

    i = np.arange(len(f.flatten())).reshape(f.shape)
    ones = np.ones_like(i).astype(float)

    # inner: laplacian
    dc = -4. * ones.copy()
    dxp = ones.copy()
    dxm = ones.copy()
    dyp = ones.copy()
    dym = ones.copy()

    ic = i
    ixp = np.roll(i, -1, axis=1)
    ixm = np.roll(i, 1, axis=1)
    iyp = np.roll(i, -1, axis=0)
    iym = np.roll(i, 1, axis=0)

    # boundary conditions: derivatives
    # laplacian = (xp - c) - (c - xm) + (yp - c) - (c - ym)
    # subtract gradients at the sides
    #dxp[:,-1] -= 1. ; dc[:,-1] += 1.
    #dxm[:, 0] -= 1. ; dc[:, 0] += 1.
    dyp[-1,:] -= 1. ; dc[-1,:] += 1.
    dym[ 0,:] -= 1. ; dc[ 0,:] += 1.

    # add given values
    #f[:,-1] += (v[:,-1] + v[:,-2]) * 0.5
    #f[:, 0] -= (v[:, 0] + v[:, 1]) * 0.5
    f[-1,:] += -(u[-1,:] + u[-2,:]) * 0.5
    f[ 0,:] -= -(u[ 0,:] + u[ 1,:]) * 0.5

    # combine
    ld = (dc, dxp, dxm, dyp, dym)
    data = np.stack(ld, axis=-1)
    li = (ic, ixp, ixm, iyp, iym)
    indices = np.stack(li, axis=-1)
    rowsize = len(li) * np.ones_like(ic)

    # flatten
    data = data.flatten()
    indices = indices.flatten()
    rowsize = rowsize.flatten()
    indptr = np.concatenate(([0], np.cumsum(rowsize)))



eps = 0.55
nx = 32

hx = 1. / nx
pi = np.pi
cos = np.cos

def eta(x):
    e = eps
    c = cos
    return (e*c(2*pi*x) + 0.5*e**2*c(4*pi*x) + 3./8.*e**3*c(6*pi*x)) / (2*pi)


x = hx*0.5 + np.linspace(-0.5, 0.5, nx, endpoint=False)

plt.plot(x, eta(x))
ax = plt.gca()
ax.set_aspect('equal')
plt.savefig('a.pdf')
