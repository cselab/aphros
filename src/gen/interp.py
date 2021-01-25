#!/usr/bin/env python3

from sympy import *
from sympy.codegen.rewriting import ReplaceOptim

x = Symbol('x')
dw = Symbol('dw')
uw = Symbol('uw')
dq = Symbol('dq')
uq = Symbol('uq')
uf = Symbol('uf')
data = [(0, uf), (dw, uw), (dq, uq)]
#data = [(0, uf), (dw, u)]
p = polys.polyfuncs.interpolate(data, x)
p = -p.diff(x)
p = factor(p.subs(x, 0))

x = Symbol('x')
sqr = Function('sqr')
sqr_custom = ReplaceOptim(lambda p: p.is_Pow and p.exp == 2,
                          lambda p: sqr(p.base))
p = sqr_custom(p)

c = printing.ccode(p, user_functions={'sqr': 'sqr'})
print(c)
