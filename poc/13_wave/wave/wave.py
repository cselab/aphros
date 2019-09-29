#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import numpy as np
import scipy
import scipy.interpolate
import scipy.sparse as sp
import scipy.sparse.linalg
import matplotlib.pyplot as plt

def report(re, n="residual"):
    re = re.flatten()
    nr = np.linalg.norm
    l = len(re)
    print("{:}: L1: {:}, L2: {:}, Linf: {:}".format(
        n, nr(re, ord=1) / l, nr(re, ord=2) / (l**0.5), nr(re, ord=np.inf)))



# xx,yy: functions of uu,vv
def CalcJ(xx, yy):
    Jxu = (np.roll(xx, -1, axis=1) - np.roll(xx, 1, axis=1)) * 0.5;
    Jxv = (np.roll(xx, -1, axis=0) - np.roll(xx, 1, axis=0)) * 0.5;
    Jyu = (np.roll(yy, -1, axis=1) - np.roll(yy, 1, axis=1)) * 0.5;
    Jyv = (np.roll(yy, -1, axis=0) - np.roll(yy, 1, axis=0)) * 0.5;
    Jxu[:,0] = Jxu[:,1]
    Jxu[:,-1] = Jxu[:,-2]
    Jyu[:,0] = Jyu[:,1]
    Jyu[:,-1] = Jyu[:,-2]
    Jxv[0,:] = Jxv[1,:]
    Jxv[-1,:] = Jxv[-2,:]
    Jyv[0,:] = Jyv[1,:]
    Jyv[-1,:] = Jyv[-2,:]
    return np.array([[Jxu, Jxv], [Jyu, Jyv]])

# xx,yy: functions of uu,vv
def CalcI(J):
    return np.linalg.inv(J.T).T

def CalcDp(s):
    su = (np.roll(s, -1, axis=1) - s)
    sv = (np.roll(s, -1, axis=0) - s)
    sv[-1,:] *= 0
    return su,sv

def CalcDm(s):
    su = (s - np.roll(s, 1, axis=1))
    sv = (s - np.roll(s, 1, axis=0))
    sv[0,:] *= 0
    return su,sv

def CalcDD(s):
    su,sv = CalcDp(s)
    suu,suv = CalcDm(su)
    svu,suv = CalcDm(sv)
    return su,(suv+svu)*0.5,sv

def Lapl(suu, suv, svv, I):
    return (suu*I[0, 0] + suv*I[1, 0])*I[0, 0] + \
           (suu*I[0, 1] + suv*I[1, 1])*I[0, 1] + \
           (suv*I[0, 0] + svv*I[1, 0])*I[1, 0] + \
           (suv*I[0, 1] + svv*I[1, 1])*I[1, 1]

# solve Posson equation for
# laplace p = -vort(u,v)
# stream function p (psi) defined via u = dp/dy, v = -dp/dx
def pois(f, itmax=10000, tol=1e-5):
    # indexing: u[y,x]
    s = np.zeros_like(f)
    it = 0
    r = tol + 1.
    itm = 0
    while r > tol and it < itmax:
        sm = s
        k = 1.0
        d = (np.roll(s, -1, axis=0) + np.roll(s, 1, axis=0) +
                np.roll(s, -1, axis=1) + np.roll(s, 1, axis=1)
                - 4. * s - f) / 4.
        s = sm + k * d
        # boundary conditions: derivatives
        gy0 = 1
        gy1 = 0
        s[0,:] = s[1,:] - gy0
        s[-1,:] = s[-2,:] + gy1
        # center at zero
        s -= s.mean()
        r = ((s - sm)**2).mean()**0.5
        # report
        if it >= itm + 100:
            print("i={:06d}, r={:.5e}".format(it, r))
            itm = it
        it += 1
    return s

def eta(x):
    e = eps
    c = cos
    return (e*c(2*pi*x) + 0.5*e**2*c(4*pi*x) + 3./8.*e**3*c(6*pi*x)) / (2*pi)

eps = 0.55
nx = 32

# range
# x=[-0.5,0.5], y=[0,1]
# x=i*hx-0.5, y=j*hx

yc = 0.5
U = 1      # wave velocity
hx = 1. / nx
pi = np.pi
cos = np.cos

x = hx*0.5 + np.linspace(-0.5, 0.5, nx, endpoint=False)
u = np.linspace(0, 1, nx, endpoint=False)
v = np.linspace(0, 1, nx, endpoint=False)

ywave = yc + eta(x)

uu,vv = np.meshgrid(u,v)

xx = uu - 0.5
yy = vv * (yc + eta(xx))

f = vv*1

#s = pois(f)
J = CalcJ(xx, yy)
I = CalcI(J)
s = Lapl(*CalcDD(f), I)
print(s.min(), s.max())

fig = plt.figure()

if 0:
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(xx, yy, s, linewidth=0, antialiased=False)
else:
    ax = fig.gca()
    plt.plot(x, ywave)
    plt.scatter(xx, yy, c=s, s=12)
    ax.set_aspect('equal')
    plt.xlim(-0.5, 0.5)
    plt.ylim(0, 1)


plt.savefig('a.pdf')
