#!/usr/bin/env python3

from sympy import *
from sympy.printing import ccode

r = S("r") # particle radius
rho_p = S("rho_p") # particle density
rho = S("rho") # fluid density
mu = S("mu") # fluid viscosity
dt = S("dt") # time step
u = S("u") # fluid velocity
v = S("v") # particle velocity
g = S("g") # gravity
vn = S("vn") # particle velocity at the next time step

vol = Rational(4, 3) * pi * r ** 3

lhs = rho_p * vol * (vn - v) / dt
rhs = 6 * pi * mu * r * (u - vn) + (rho_p - rho) * g * vol

expr = solve(lhs - rhs, vn)[0]
expr = collect(expr, g)
pprint([lhs, rhs])
pprint(expr)

user_functions = {'sqr': 'sqr'}
user_functions = {
    'Pow': [(lambda b, e: e == 2, lambda b, e: 'sqr(%s)' % b),
            (lambda b, e: e != 2, 'pow')]
}
code = ccode(expr, assign_to="v", user_functions=user_functions)
print(code)
