#!/usr/bin/env python

from sympy import *
from sympy.printing import ccode
import textwrap


u,v = S("u,v")
X = Function('X')(u,v)
Y = Function('Y')(u,v)
s = Function('s')

J_ = Matrix([X, Y]).jacobian(Matrix([u, v]))
J = MatrixSymbol('J', 2, 2)
I_ = J.inv()
I = MatrixSymbol('I', 2, 2)

def DX(s):
    return diff(s,u) * I[0,0] + diff(s,v) * I[1,0]

def DY(s):
    return diff(s,u) * I[0,1] + diff(s,v) * I[1,1]

def LAPL(s):
    return DX(DX(s)) + DY(DY(s))

s =  Function('s')(u,v)

#pprint(LAPL(s))

suu,suv,svu,svv = S('suu,suv,svu,svv')
suu_ = diff(s,u,u)
suv_ = diff(s,u,v)
svu_ = diff(s,v,u)
svv_ = diff(s,v,v)
sd_ = [suu_, suv_, svu_, svv_]
sd = [suu, suv, svu, svv]

pprint(LAPL(s).subs(zip(sd_, sd)))
