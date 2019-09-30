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
def CalcJp(xx, yy):
    Jxu = np.roll(xx, -1, axis=1) - xx
    Jxv = np.roll(xx, -1, axis=0) - xx
    Jyu = np.roll(yy, -1, axis=1) - yy
    Jyv = np.roll(yy, -1, axis=0) - yy
    Jxu[:,-1] = Jxu[:,-2]
    Jyu[:,-1] = Jyu[:,-2]
    Jxv[-1,:] = Jxv[-2,:]
    Jyv[-1,:] = Jyv[-2,:]
    return np.array([[Jxu, Jyu], [Jxv, Jyv]])

def CalcJm(xx, yy):
    Jxu = xx - np.roll(xx, 1, axis=1)
    Jxv = xx - np.roll(xx, 1, axis=0)
    Jyu = yy - np.roll(yy, 1, axis=1)
    Jyv = yy - np.roll(yy, 1, axis=0)
    Jxu[:,0] = Jxu[:,1]
    Jyu[:,0] = Jyu[:,1]
    Jxv[0,:] = Jxv[1,:]
    Jyv[0,:] = Jyv[1,:]
    return np.array([[Jxu, Jyu], [Jxv, Jyv]])

# xx,yy: functions of uu,vv
def CalcI(J):
    return np.linalg.inv(J.T).T

def CalcDp(s, Ip):
    I=Ip
    su = (np.roll(s, -1, axis=1) - s)
    sv = (np.roll(s, -1, axis=0) - s)
    su[:,-1] = su[:,-2]
    sv[-1,:] = sv[-2,:]
    return I[0,0]*su+I[0,1]*sv, I[1,0]*su+I[1,1]*sv

def CalcDm(s, Im):
    I=Im
    su = (s - np.roll(s, 1, axis=1))
    sv = (s - np.roll(s, 1, axis=0))
    su[:,0] = su[:,1]
    sv[0,:] = sv[1,:]
    return I[0,0]*su+I[0,1]*sv, I[1,0]*su+I[1,1]*sv

def CalcDc(s, Im, Ip):
    mu,mv = CalcDm(s, Im)
    pu,pv = CalcDp(s, Ip)
    return (mu+pu)*0.5, (mv+pv)*0.5


def CalcQDc(s):
    I=Im
    su = (np.roll(s, -1, axis=1) - np.roll(s, 1, axis=1)) * 0.5
    sv = (np.roll(s, -1, axis=0) - np.roll(s, 1, axis=0)) * 0.5
    su[:,0] = su[:,1]
    sv[0,:] = sv[1,:]
    su[:,-1] = su[:,-2]
    sv[-1,:] = sv[-2,:]
    return su, sv

def CalcDD(s, Im, Ip):
    su,sv = CalcDp(s, Ip)
    suu,suv = CalcDm(su, Im)
    svu,svv = CalcDm(sv, Im)
    return suu,(suv+svu)*0.5,svv

def Lapl(s, Im, Ip):
    suu,suv,svv = CalcDD(s, Im, Ip)
    return suu + svv

# solve Posson equation
# U: ds/dx=U*nx
def pois(f, Im, Ip, Jm, Jp, U, itmax=10000, tol=1e-5, relax=0.3):
    # indexing: u[y,x]
    s = np.zeros_like(f)
    it = 0
    r = tol + 1.
    itm = 0
    while r > tol and it < itmax:
        sm = s
        k = relax
        Imd = Im[0,0]*Im[1,1] - Im[1,0]*Im[0,1]
        Ipd = Ip[0,0]*Ip[1,1] - Ip[1,0]*Ip[0,1]
        Id = (Imd + Ipd) * 0.5
        d = (Lapl(s, Im, Ip) - f) / Id / 4
        s = sm + k * d
        # derivative in u
        su = CalcQDc(s)[0][-1,:]

        I = (Im + Ip) * 0.5
        J = (Jm + Jp) * 0.5
        # tangent to boundary
        tx = J[0,0][-1,:]
        ty = J[0,1][-1,:]
        # outer normal to boundary
        nx = -ty
        ny = tx
        nn = (nx**2 + ny**2)**0.5
        nx /= nn
        ny /= nn

        # derivative in v
        sv = (U*nx - (I[0,0][-1,:]*nx+I[0,1][-1,:]*ny)*su) /  \
             (I[1,0][-1,:]*nx+I[1,1][-1,:]*ny)

        s[0,:] = s[1,:]
        s[-1,:] = s[-2,:] + sv
        # center at zero
        s -= s.max()
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

if 0:
    a = pi*0.125
    sa = np.sin(a)
    ca = np.cos(a)
    xx = (uu*ca + -vv*sa)*0.5
    yy = (uu*sa + vv*ca)*0.5

f = xx*0

Jp = CalcJp(xx, yy)
Jm = CalcJm(xx, yy)
Im = CalcI(Jm)
Ip = CalcI(Jp)

s = pois(f, Im, Ip, Jm, Jp, U)

s = CalcDc(s, Im, Ip)[1]

#fig,ax = plt.figure()
fig,(ax1,ax2) = plt.subplots(1,2)

if 0:
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(xx, yy, s, linewidth=0, antialiased=False)
else:
    ax = fig.gca()
    #ax1.plot(x, ywave)
    ax1.scatter(xx, yy, c=f, s=20)
    ax1.set_xlim(-0.5, 0.5)
    ax1.set_ylim(0, 1)
    ax1.set_aspect('equal')

    ax2.scatter(xx, yy, c=s, s=20)
    ax2.set_xlim(-0.5, 0.5)
    ax2.set_ylim(0, 1)
    ax2.set_aspect('equal')


plt.savefig('a.pdf')
