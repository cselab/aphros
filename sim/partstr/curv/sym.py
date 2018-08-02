#!/usr/bin/env python

import sympy as sp
import sympy.vector
import numpy as np

def Curv():
    u = sp.Symbol('u')
    v = sp.Symbol('v')

    a  = 0.54353
    b  = 0.121435
    c  = -0.561365

    x = sp.cos(v) * sp.cos(u)
    y = sp.cos(v) * sp.sin(u)
    q = sp.sin(v)
    z = a * q + b * q ** 3 + c * q ** 5

    r = sp.Matrix([x, y, z])
    ru = sp.diff(r, u)
    rv = sp.diff(r, v)

    E = ru.dot(ru)
    F = ru.dot(rv)
    G = rv.dot(rv)

    n = ru.cross(rv)
    n /= sp.sqrt(n.dot(n))

    ruu = sp.diff(ru, u)
    ruv = sp.diff(ru, v)
    rvv = sp.diff(rv, v)

    e = ruu.dot(n)
    f = ruv.dot(n)
    g = rvv.dot(n)

    H = (e * g - f * f) / (E * G - F * F)

    H = sp.simplify(H)

    fH = sp.lambdify([u, v], H)
    fx = sp.lambdify([u, v], x)
    fy = sp.lambdify([u, v], y)
    fz = sp.lambdify([u, v], z)
    return fx,fy,fz,fH


fx,fy,fz,fH = Curv()

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

u = np.linspace(-np.pi, np.pi, 50)
v = np.linspace(-np.pi * 0.5, np.pi * 0.5, 50)
u, v = np.meshgrid(u, v)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(u, v, fH(u, v))
ax.plot_surface(u, v, fz(u, v))
#ax.plot_surface(fx(u,v), fy(u,v), fz(u,v))
plt.show()


