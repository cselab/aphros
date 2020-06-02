#!/usr/bin/env python3

from sympy import *
from itertools import product
from sympy.printing import ccode
import textwrap

dh = S.One / 2

#u = Function('u')
u = lambda x,y,z : S("u({:},{:},{:})".format(x,y,z))

# interpolation to node
# u: function defined on cells
# x,y,z: node index
def tonode(u, x, y, z):
    r = S('0')
    for sx,sy,sz in product([-1,1], repeat=3):
        r = r + u(x + sx * dh, y + sy * dh, z + sz * dh)
    return r / 8


# average over nodes
# u: function defined on nodes
# x,y,z: cell index
def average(u, x, y, z):
    r = x*0
    w = 0
    for sx,sy,sz in product([-1,1], repeat=3):
        r += u(x + sx * dh,
               y + sy * dh,
               z + sz * dh)
        w += 1
    return r / w

e = average(lambda x,y,z: tonode(u, x, y, z), 0, 0, 0)
e = simplify(e)

s = []
for sx,sy,sz in product([-1,0,1], repeat=3):
    t = 'u({:},{:},{:})'.format(sz,sy,sx)
    s += [float(e.coeff(t))]


o = "std::array<Scal, SN> a = {\nST};"
o = o.replace('SN', str(len(s)))
o = o.replace('ST', ', '.join(map(str, s)))
o = textwrap.fill(o, 60, subsequent_indent='    ')
print(o)

