#!/usr/bin/env python

'''
Estimation of curvature using connected particles.

Petr Karnakov (kpetr@ethz.ch)
2019-06-01
ETH Zurich
'''

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Common notation:
# xx: particle string, rows (x,y), shape (N,2)
# ff: forces, rows (x,y),  shape (N,2)
# ll: interface lines, rows (x0,y0,x1,y1), shape (N, 4)
# p: position of central particle
# ph: phi, orientation angle
# th: theta, bending angle

# Global:
N = 9                # particles per string
c = N // 2           # index of central particle
hp = 4. / (N - 1)    # distance between particles

# Point from coordinates.
def P(x0, x1):
    return np.array((x0, x1))

# Scalar product.
def Dot(x, y):
    return (x * y).sum()

def Norm(x):
    return Dot(x, x) ** 0.5

def Dist(x, y):
    return Norm(x - y)

# Unit vector at angle ph.
def E(ph):
    return P(np.cos(ph), np.sin(ph))

# Positions of particles.
def X(p, ph, th):
    global N, c, hp
    xx = np.zeros((N,2))
    xx[c] = p
    for j in range(0, c):
        jp = j + 0.5
        xx[c + j + 1] = xx[c + j] + hp * E(ph + jp * th)
        xx[c - j - 1] = xx[c - j] - hp * E(ph - jp * th)
    return xx

# Derivative of positions by phi.
def DxDph(p, ph, th):
    global N, c, hp
    xx = np.zeros((N,2))
    for j in range(0, c):
        jp = j + 0.5
        xx[c + j + 1] = xx[c + j] + hp * E(ph + jp * th + np.pi * 0.5)
        xx[c - j - 1] = xx[c - j] - hp * E(ph - jp * th + np.pi * 0.5)
    return xx

# Dereivative of positions by theta.
def DxDth(p, ph, th):
    global N, c, hp
    xx = np.zeros((N,2))
    for j in range(0, c):
        jp = j + 0.5
        xx[c + j + 1] = xx[c + j] + hp * jp * E(ph + jp * th + np.pi * 0.5)
        xx[c - j - 1] = xx[c - j] + hp * jp * E(ph - jp * th + np.pi * 0.5)
    return xx

# Nearest to x point on line [a,b].
# a,b: line endpoints, shape=(2)
# x: point, shape=(2)
def Nearest(a, b, x):
    q = b - a
    k = Dot(q, x - a) / Dot(q, q)
    k = np.clip(k, 0., 1.)
    return a + q * k;

# Displacement from line to circular arc.
# k: curvature
# w: distance from endpoint to center
# d: distance from point to center
def Delta(k, w, d):
    k = min(k, 1. / w)      # limit radius of curvature by w
    d = min(d, w)
    t1 = (1. - (k * w) ** 2) ** 0.5
    t2 = (1. - (k * d) ** 2) ** 0.5
    return k * (w ** 2 - d ** 2) / (t1 + t2);

# Point on circular arc.
# k: curvature of circular arc
# a, b: endpoints
# y: point on line
def Circ(k, a, b, y):
    q = b - a
    n = P(q[1], -q[0]) / Norm(q)   # outer unit normal
    c = (a + b) * 0.5
    d = Dist(c, y)
    w = Dist(c, a)
    return y + Delta(k, w, d) * n;

# Point on circular arc at nearest line.
# ll: list of lines (x0,y0,x1,y1)
# x: point
# k: curvature of circular arc
def NearestCirc(ll, x, k):
    yy = [Nearest(l[0:2], l[2:4], x) for l in ll]
    i = np.argmin([Dist(x, y) for y in yy])
    y = yy[i]
    a = ll[i,0:2]
    b = ll[i,2:4]
    xl = Circ(k, a, b, y)
    return xl

# Force.
# ll: lines
# xx: particle positions
# k: curvature
# eta: relaxation
def F(ll, xx, k, eta):
    ff = [eta * (NearestCirc(ll, x, k) - x) for x in xx]
    return ff

# Parameters after one iteration.
def Iter(p, ph, th, ff):
    p = np.copy(p)
    ff = np.copy(ff)

    xx = X(p, ph, th)

    # correct p
    p += ff[c]
    dx = X(p, ph, th) - xx
    xx += dx
    ff -= dx

    # correct phi
    d = DxDph(p, ph, th)
    ph += Dot(ff, d) / Dot(d, d)
    dx = X(p, ph, th) - xx
    xx += dx
    ff -= dx

    # correct theta
    d = DxDth(p, ph, th)
    th += Dot(ff, d) / Dot(d, d)
    dx = X(p, ph, th) - xx
    xx += dx
    ff -= dx

    return (p, ph, th)

# Curvature from theta.
def Curv(th):
    return 2 ** 0.5 / hp * np.sin(th) / (1. + np.cos(th)) ** 0.5

# p,ph,th: initial
# ll: lines
# eta: relaxation
# eps: tolerance
# mmax: maximum number of iterations
# h: cell size
def Solve(p, ph, th, ll, eta=0.5, eps=1e-5, mmax=20, h=1.):
    xx = X(p, ph, th)
    for m in range(mmax):
        (p, ph, th) = Iter(p, ph, th, F(ll, xx, Curv(th), eta))
        dx = X(p, ph, th) - xx
        r = abs(dx).max() / (eta * h)
        if r < eps:
            break
        xx += dx
    return (p, ph, th)

def PlotLines(ax, ll):
    for i in range(ll.shape[0]):
        x0 = ll[i][np.r_[0,2]]
        x1 = ll[i][np.r_[1,3]]
        ax.plot(x0, x1, c="black")

# Lines with endpoints on unit circle.
# n: number of lines
def Lines(n):
    np.random.seed(1)
    q = 2. * np.pi / n     # increment of angle
    dq = q * 0.5           # noise magnitude
    rr = np.random.rand(n) * dq
    rr = np.append(rr, rr[0])
    ll = [np.hstack((
            E(q * i + rr[i]),
            E(q * (i + 1) + rr[i + 1]))) for i in range(n)]
    ll = np.array(ll)
    return ll

# c0, c1: parameters, (p, ph, th)
# o: path to output
def Plot(c0, c1, o):
    fig,ax = plt.subplots()

    ax.set_aspect("equal")
    e = 3
    ax.set_xlim(-e, e)
    ax.set_ylim(-e, e)

    PlotLines(ax, ll)

    xx = X(*c0)
    plt.scatter(*xx.T, c="red")
    xx = X(*c1)
    plt.scatter(*xx.T)

    plt.savefig("a.pdf")


ll = Lines(7)

# initial
a = ll[0][0:2]
b = ll[0][2:4]
p = (a + b) * 0.5
ph = np.arctan2((b - a)[1], (b - a)[0])
th = 0.
c0 = (p, ph, th)

# final
c1 = Solve(c0[0], c0[1], c0[2], ll, eps=1e-5, mmax=20)

# curvature
print(Curv(c1[2]))

# plot initial and final
Plot(c0, c1, "a.pdf")
