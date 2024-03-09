#!/usr/bin/env python3

import pint
import math


def pprint(expr, desc=None):
    res = eval(expr, globals()).to_base_units()
    print("{} = {:.6g}".format(
        '    ' + expr if not desc else desc + ':\n    ' + expr, res))


ureg = pint.UnitRegistry()

R = 8.314 * ureg.J / (ureg.K * ureg.mol)  # gas constant
T = 293 * ureg.K  # temperature
p = 10**5 * ureg.Pa  # pressure
mol_to_c = R * T / p  # molar concentration to dimensionless concentration
k_H2 = 7.7e-6 * ureg.mol / (ureg.m**3 * ureg.Pa)  # Henry coefficient
lambda_H2 = k_H2 * R * T  # Ostwald coefficient

I = 10e-6 * ureg.A  # current
ri = 400e-6 * ureg.m  # cathode inner radius
re = 410e-6 * ureg.m  # cathode outer radius
A = math.pi * (re**2 - ri**2)  # cathode area
n = 2  # stoichiometric constant for H2
F = 96485 * ureg.C / ureg.mol  # Faraday constant
D = 5e-9 * ureg.m**2 / ureg.s  # diffusion coefficient
t_D = re**2 / D  # diffusion time
sigma = 72e-3 * ureg.N / ureg.m
L = re
mu = 1e-3 * ureg.Pa * ureg.s
rho = 1000 * ureg.kg / ureg.m**3
nu = mu / rho
#E = 2.1 * ureg.V  # typical cell potential
#Rcell = E / I  # cell resistance
#vel_D = re / t_D  # diffusion velocity
#g = 9.8 * ureg.m / ureg.s**2

charge_to_mol = 1 / (n * F)  # charge to amount of substance of dissolved gas
charge_to_c = charge_to_mol * mol_to_c

Oh = mu / (rho * L * sigma)**0.5
Sc = nu / D
#Ca = mu * vel_D / sigma
#Re = L * vel_D / nu
#Eo = rho * g * L**2 / sigma
#Ga = g * L**3 / nu**2
#g_tilde = g / L * t_D**2

molflux = (I / A) / (n * F)
molrate = molflux * A  # molar rate of dissolved gas

cflux = molflux * mol_to_c
cflux_tilde = cflux / (L / t_D)

crate = molrate * mol_to_c
crate_tilde = crate / (L**3 / t_D)

pprint("mol_to_c", "molar concentration to dimensionless concentration")
pprint("charge_to_mol", "charge to amount of substance of dissolved gas")
pprint("charge_to_c", "charge to volume of dissolved gas")
pprint("A", "current")
pprint("I/A", "current density on cathode")
pprint("molflux", "molar flux of dissolved gas")
pprint("molrate", "molar rate of dissolved gas")
pprint("cflux", "concentration flux of dissolved gas")
pprint("crate", "concentration rate of dissolved gas")
pprint("t_D", "diffusion time")
#pprint("vel_D", "diffusion velocity")
#pprint("E", "cell potential")
#pprint("Rcell", "cell resistance")
pprint("sigma", "surface tension")
pprint("k_H2", "Henry coefficient")
pprint("lambda_H2", "Ostwald coefficient")
#pprint("g", "gravity")

print("\nDimensionless:")
pprint("cflux_tilde", "molar flux of dissolved gas")
pprint("crate_tilde", "molar rate of dissolved gas")
#pprint("g_tilde", "gravity")
#pprint("Ca", "Capillary number")
#pprint("Re", "Reynolds number")
#pprint("Eo", "Eötvös number")
#pprint("Ga", "Galilei number")
pprint("Oh", "Ohnesorge number")
pprint("Sc", "Schmidt number")
