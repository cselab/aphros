#!/usr/bin/env python3

import pint
import math


def pprint(expr, desc=None, compact=False):
    res = eval(expr, globals()).to_base_units()
    if compact:
        res = res.to_compact()
    print("{} = {:.6g}".format(
        '    ' + expr if not desc else desc + ':\n    ' + expr, res))


ureg = pint.UnitRegistry()

R = 8.314 * ureg.J / (ureg.K * ureg.mol)  # gas constant
T = 293 * ureg.K  # temperature
p = 10**5 * ureg.Pa  # pressure
mol_to_c = R * T / p  # molar concentration to dimensionless concentration
k_H2 = 7.7e-6 * ureg.mol / (ureg.m**3 * ureg.Pa)  # Henry coefficient
lambda_H2 = k_H2 * R * T  # Ostwald coefficient

Q = 100 * ureg.ml / ureg.min  # flow rate
depth = 1.5 * ureg.mm
inlet_width = 1.5 * ureg.mm
inlet_area = depth * inlet_width
inlet_velocity = Q / inlet_area
device_width = 8.8 * ureg.mm

J = 700 * ureg.mA / ureg.cm**2  # current density
n = 2  # stoichiometric constant for H2
F = 96485 * ureg.C / ureg.mol  # Faraday constant
D_H2 = 6e-9 * ureg.m**2 / ureg.s  # diffusion coefficient for H2
D_O2 = 2.5e-9 * ureg.m**2 / ureg.s  # diffusion coefficient for O2
t_unit = inlet_width / inlet_velocity  # unit time
L_unit = device_width  # unit length
sigma = 72e-3 * ureg.N / ureg.m
mu = 1e-3 * ureg.Pa * ureg.s
rho = 1000 * ureg.kg / ureg.m**3
nu = mu / rho  # kinematic viscosity
E = 20 * ureg.V  # typical cell potential
Rcell = E / J  # cell resistivity

charge_to_mol = 1 / (n * F)  # charge to amount of substance of dissolved gas
charge_to_c = charge_to_mol * mol_to_c

Sc_H2 = nu / D_H2
Sc_O2 = nu / D_O2
Ca = mu * inlet_velocity / sigma
Re = inlet_width * inlet_velocity / nu

molflux = J / (n * F)

cflux = molflux * mol_to_c
cflux_tilde = cflux / (L_unit / t_unit)

inlet_width_tilde = inlet_width / L_unit

pprint("mol_to_c", "molar concentration to dimensionless concentration")
pprint("charge_to_mol", "charge to amount of substance of dissolved gas")
pprint("charge_to_c", "charge to volume of dissolved gas")
pprint("J", "current density on cathode")
pprint("molflux", "molar flux of dissolved gas")
pprint("cflux", "concentration flux of dissolved gas")
pprint("t_unit", "time unit")
pprint("E", "cell potential")
pprint("Rcell", "cell resistivity")
pprint("sigma", "surface tension")
pprint("k_H2", "Henry coefficient")
pprint("lambda_H2", "Ostwald coefficient")
pprint("Q", "flowrate")
pprint("inlet_velocity", "inlet mean velocity")

print("\nDimensionless:")
pprint("cflux_tilde", "molar flux of dissolved gas")
pprint("inlet_width_tilde", "inlet width")
pprint("Ca", "Capillary number")
pprint("Re", "Reynolds number")
pprint("Sc_H2", "Schmidt number H2")
pprint("Sc_O2", "Schmidt number O2")
