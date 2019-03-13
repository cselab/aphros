#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np

ll = np.loadtxt(sys.argv[1])
n  = int(sys.argv[2])

nh = n // 2
eta = 0.5 # relaxation
hp = 4. / (n - 1)
h = 1.
sc = 1.
pi = np.pi
itermax = 1000

# xx: particle string, rows (x,y), shape (n,2)
# ff: forces, rows (x,y),  shape (n,2)
# ll: interface lines, rows (x0,y0,x1,y1), shape (n, 4)
# p: origin
# a: alpha
# t: theta

def P(x0, x1):
    return np.array((x0, x1))

def Dot(x, y):
    return (x * y).sum()

# x: point, shape=(2)
def Norm(x):
    return Dot(x, x) ** 0.5

def Dist(x, y):
    return Norm(x - y)

def E(a):
    return P(np.cos(a), np.sin(a))

def X(p, a, t):
    xx = np.zeros((n,2))
    xx[nh] = p
    for j in range(0, nh):
        jp = j + 0.5
        xx[nh + j + 1] = xx[nh + j] + hp * E(a + t * jp)
        xx[nh - j - 1] = xx[nh - j] - hp * E(a - t * jp)
    return xx

def Dxda(p, a, t):
    xx = np.zeros((n,2))
    for j in range(0, nh):
        jp = j + 0.5
        xx[nh + j + 1] = xx[nh + j] + hp * E(a + t * jp + pi * 0.5)
        xx[nh - j - 1] = xx[nh - j] - hp * E(a - t * jp + pi * 0.5)
    return xx

def Dxdt(p, a, t):
    xx = np.zeros((n,2))
    for j in range(0, nh):
        jp = j + 0.5
        xx[nh + j + 1] = xx[nh + j] + hp * jp * E(a + t * jp + pi * 0.5)
        xx[nh - j - 1] = xx[nh - j] + hp * jp * E(a - t * jp + pi * 0.5)
    return xx


# x0,x1: line ends, shape=(2)
# x: point, shape=(2)
def Nearest(x0, x1, x):
    d = x1 - x0
    k = Dot(d, x - x0) / Dot(d, d)
    k = np.clip(k, 0., 1.)
    return x0 + d * k;

# k: curvature
# l: distance from end to center
# d: distance from point to center
def SegCirc(k, l, d):
    t1 = (1. - (k * l) ** 2) ** 0.5
    t2 = (1. - (k * d) ** 2) ** 0.5
    return k * (l ** 2 - d ** 2) / (t1 + t2);

def ShSegCirc(k, x0, x1, x):
    # outer unit normal
    ld = x1 - x0
    ln = P(ld[1], -ld[0]) / Norm(ld)
    # center
    xc = (x0 + x1) * 0.5;
    # distance from center
    dc = Dist(xc, x)
    # max distance from center
    mdc = Dist(xc, x0)
    # shift from line to circle
    s = sc * SegCirc(k, mdc, dc);
    return x + ln * s;

def NearestToLines(ll, x, k):
    yy = [Nearest(l[0:2], l[2:4], x) for l in ll]
    dd = np.array([Dist(x, y) for y in yy])
    i = dd.argmin()
    y = yy[i]
    x0 = ll[i,0:2]
    x1 = ll[i,2:4]
    y = ShSegCirc(k, x0, x1, y)
    return y

def Force(ll, xx, k):
    ff = [eta * (NearestToLines(ll, x, k) - x) for x in xx]
    return ff

def Step(p, a, t, ff, cf=True):
    xx = X(p, a, t)

    # correct p
    p += ff[nh]
    dx = X(p, a, t) - xx
    if cf: ff -= dx
    xx += dx

    # correct alpha
    d = Dxda(p, a, t)
    a += Dot(ff, d) / Dot(d, d)
    dx = X(p, a, t) - xx
    if cf: ff -= dx
    xx += dx

    # correct theta
    d = Dxdt(p, a, t)
    t += Dot(ff, d) / Dot(d, d)
    dx = X(p, a, t) - xx
    if cf: ff -= dx
    xx += dx

    return (p, a, t)

def Curv(t):
    return 2 ** 0.5 / hp * np.sin(t) / (1. + np.cos(t)) ** 0.5

def Solve(p, a, t, ll, eps):
    xx = X(p, a, t)
    for i in range(itermax):
        (p, a, t) = Step(p, a, t, Force(ll, xx, Curv(t)), False)
        dx = X(p, a, t) - xx
        r = abs(dx).max() / (eta * h)
        if r < eps:
            break
        xx += dx
    return (p, a, t)

def PlotInterface(ax, ll):
    for i in range(ll.shape[0]):
        x0 = ll[i][np.r_[0,2]]
        x1 = ll[i][np.r_[1,3]]
        ax.plot(x0, x1, c="black")

p = P(0.,ll[:,np.r_[1,3]].max())
a = pi
t = 0.
p,a,t = Solve(p, a, t, ll, eps=1e-5)
xx = X(p, a, t)

for x in xx:
    print(x[0], x[1])
