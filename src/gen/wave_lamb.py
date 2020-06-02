#!/usr/bin/env python3

# Capillary effects on wave breaking
# Luc Deike, Stephane Popinet and W. Kendall Melville
# https://doi.org/10.1017/jfm.2015.103

from sympy import *
from sympy.printing import ccode
import textwrap


x,y = S("x,y")
a = S('a')
k = S('k')
h = S('h')   # bottom
g = S('g')   # gravity


eps_ = a * k

kx = k * x

chi_ = 1 / tanh(k*h)

chi = S('chi')
eps = S('eps')

# interface heigth
eta = (
    a * cos(kx) +
    eps * a * chi * (3*chi**2 - 1) * cos(2*kx) / 4 +
    eps**2*a* (
      -3*(chi**4 - 3*chi**2 + 3) * cos(kx)/8 +
      + 3*(8*chi**6+(chi**2-1)**2)*cos(3*kx)/64
    )
  )

# velocity potential

omega_ = sqrt(g*k*tanh(k*h)*(1+eps**2*(9*(chi**2-1)**2/8+chi**2)))
omega = S('omega')
phi = (
    a*g/omega*cosh(k*(y+h))/cosh(k*h)*sin(kx)+
    eps*3*a*g/omega*(chi**2-1)**2/(8*chi)*cosh(2*k*(y+h))/cosh(2*k*h)*sin(2*kx)+
        eps**2*a*g/(64*omega)*(chi**2-1)*(chi**2+3)*(9*chi**2-13)*
        cosh(3*k*(y+h))/cosh(3*k*h)*sin(3*kx)
  )

vx = diff(phi, x)
vy = diff(phi, y)

o = []

def A(e, name):
    global o
    t = ccode(e, assign_to=name)
    o += [textwrap.fill(t, 60)]

A(eps_, "eps")
A(chi_, "chi")
A(eta, "eta")
A(omega_, "omega")
A(vx, "vx")
A(vy, "vy")

o = '\n\n'.join(o)

print(o)
