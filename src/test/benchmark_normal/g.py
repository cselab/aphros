#!/usr/bin/env python

from sympy import *
from itertools import product
from sympy.printing import ccode
import textwrap

dx = Symbol('dx') # step
dxh = dx / 2
dy = Symbol('dy') # step
dyh = dy / 2
dz = Symbol('dz') # step
dzh = dz / 2
c = Symbol('c')  # cell index

u = Function('u')
hx = Symbol('h[0]')
hy = Symbol('h[1]')
hz = Symbol('h[2]')
hh = [hx, hy, hz]

# interpolation to node
# u: function defined on cells
# x,y,z: node index
def f(u, x, y, z):
    r = u(x,y,z)*0
    for sx,sy,sz in product([-1,1], repeat=3):
        r += u(x + sx * dxh, y + sy * dyh, z + sz * dzh)
    return r / 8


# gradient from nodes
# u: function defined on nodes
# x,y,z: cell index
# d: direction, 0..2
def g(u, x, y, z, d):
    r = x*0
    for ss in product([-1,1], repeat=3):
        r += u(x + ss[0] * dxh,
               y + ss[1] * dyh,
               z + ss[2] * dzh) * ss[d] / hh[d]
    return r / 4

for d in range(3):
    e = g(lambda x,y,z: f(u, x, y, z), 0, 0, 0, d)
    e = simplify(e)
    e = e.subs({dx: 1, dy: 1, dz: 1})
    uf = {
            'u' : [(lambda x,y,z: True,
              lambda x,y,z: "q({:},{:},{:})".format(x,y,z))]
            , 'hx' : 'h[0]'
            , 'hy' : 'h[1]'
            , 'hz' : 'h[2]'
            }
    o = ccode(e, user_functions=uf, assign_to='pn[i][{:}]'.format(d))
    o = textwrap.fill(o, 60)
    print(o)
