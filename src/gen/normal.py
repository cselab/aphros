#!/usr/bin/env python3

from sympy import *
from itertools import product
from sympy.printing import ccode
import textwrap

dx = Symbol('dx')  # step
dxh = dx / 2
dy = Symbol('dy')  # step
dyh = dy / 2
dz = Symbol('dz')  # step
dzh = dz / 2
c = Symbol('c')  # cell index

u = Function('u')
hx = Symbol('h[0]')
hy = Symbol('h[1]')
hz = Symbol('h[2]')


def Run2d():
    dim = 2
    def interp(u, x, y):
        '''
        interpolation to node
        u: function defined on cells
        x,y: node index
        '''
        r = u(x, y) * 0
        for sx, sy in product([-1, 1], repeat=dim):
            r += u(x + sx * dxh, y + sy * dyh)
        return r / 8

    def grad(u, x, y, d):
        '''
        gradient from nodes
        u: function defined on nodes
        x,y: cell index
        d: direction, 0..2
        '''
        hh = [hx, hy]
        r = x * 0
        for ss in product([-1, 1], repeat=dim):
            r += u(x + ss[0] * dxh, y + ss[1] * dyh) * ss[d] / hh[d]
        return r / 4

    for d in range(dim):
        e = grad(lambda x, y: interp(u, x, y), 0, 0, d)
        e = simplify(e)
        e = e.subs({dx: 1, dy: 1})
        uf = {
            'u':
            [(lambda x, y: True, lambda x, y: "q({:},{:})".format(x, y))],
            'hx': 'h[0]',
            'hy': 'h[1]',
        }
        o = ccode(e, user_functions=uf, assign_to='pn[i][{:}]'.format(d))
        o = textwrap.fill(o, 60)
        print(o)


def Run3d():
    dim = 3
    def interp(u, x, y, z):
        '''
        interpolation to node
        u: function defined on cells
        x,y,z: node index
        '''
        r = u(x, y, z) * 0
        for sx, sy, sz in product([-1, 1], repeat=dim):
            r += u(x + sx * dxh, y + sy * dyh, z + sz * dzh)
        return r / 8

    def grad(u, x, y, z, d):
        '''
        gradient from nodes
        u: function defined on nodes
        x,y,z: cell index
        d: direction, 0..2
        '''
        hh = [hx, hy, hz]
        r = x * 0
        for ss in product([-1, 1], repeat=dim):
            r += u(x + ss[0] * dxh, y + ss[1] * dyh,
                   z + ss[2] * dzh) * ss[d] / hh[d]
        return r / 4

    for d in range(dim):
        e = grad(lambda x, y, z: interp(u, x, y, z), 0, 0, 0, d)
        e = simplify(e)
        e = e.subs({dx: 1, dy: 1, dz: 1})
        uf = {
            'u': [(lambda x, y, z: True,
                   lambda x, y, z: "q({:},{:},{:})".format(x, y, z))],
            'hx':
            'h[0]',
            'hy':
            'h[1]',
            'hz':
            'h[2]',
        }
        o = ccode(e, user_functions=uf, assign_to='pn[i][{:}]'.format(d))
        o = textwrap.fill(o, 60)
        print(o)


print("// 2D")
Run2d()
print("// 3D")
Run3d()
