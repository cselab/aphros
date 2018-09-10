#!/usr/bin/env python

import argparse
import math

# Returns dictionary with standard parameters {k : (v,h)}
# k: name
# v: default value or type [int, float, list] (then default is None)
# h: description for help
def StdParAll():
    Re = (float, "Reynolds number")
    Ca = (float, "capillary number")
    We = (float, "Weber number")
    La = (float, "Laplace number")
    Ga = (float, "Galilei number")
    np = (1, "number of processors")
    nx = (64, "mesh size")
    brh = (4., "br/h, cells per bubble radius")
    tmax = (1., "total time")
    nfr = (100, "number of frames")
    mu = (1., "dynamic viscosity of bubble relative to carrier")
    rho0 = (1., "density of carrier fluid")
    rho = (1., "density of bubble relative to carrier")
    dim = (3, "dimension, 2 or 3")
    chsm = (1, "ch smooth steps (rhor>10: 1, rhor>100: 2)")
    gesm = (1, "ge smooth steps")
    wallx = (0, "wall in x direction, else periodic")
    wally = (0, "wall in y direction, else periodic")
    wallz = (0, "wall in z direction, else periodic")
    return locals().copy()

# ordering:
# x:    length
# vel:  velocity
# rho:  density
# mu:   dynamic viscosity
# sig:  surface tension

def GetRe(x, vel, rho, mu):
    return x * vel * rho / mu

def GetCa(vel, mu, sig):
    return vel * mu / sig

def GetWe(x, vel, rho, sig):
    return rho * vel ** 2 * x / sig

def GetLa(x, rho, mu, sig):
    return rho * sig * x / mu ** 2

def GetGa(x, rho, mu, sig):
    return rho * sig * x / mu ** 2

# Returns a selection from StdParAll().
# kk: list of names
def StdPar(l):
    d = StdParAll()
    r = dict()  # result
    for k in kk:
        r[k] = d[k]
    return r;

def IsClose(a, b):
    return abs(a - b) < 1e-10

# Returns Namespace with parameters from c overriden by args
# c: dict with parameters as in StdParAll()
# desc: description for command line help
def GetArgs(c, desc):
    p = argparse.ArgumentParser(description=desc,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    for k in c:
        v,h = c[k]
        if v[0] in [int, float, list]:
            t = v
            v = None
        else:
            t = type(v)
            v = v

        n = '-' + k

        if t == list:
            m = ('X','Y','Z')
            assert v is None or len(v) == len(m)
            t = type(v[0]) if v is not None else float
            p.add_argument(n, nargs=len(m), metavar=m, default=v, help=h, type=t)
        else:
            p.add_argument(n, nargs='1', default=v, help=h, type=t)

    a = p.parse_args()
    return a

# Returns python code defining variables from dict c
def GetDictPy(c):
    t = ""
    for k in c:
        v = c[k]
        v = '"{:}"'.format(v) if type(v) == str else str(v)
        t += "{:} = {:}\n".format(k, v)

    return t

def norm(v):
    assert len(v) == 3
    return sum([a ** 2 for a in v]) ** 0.5

# Returns dictionary with derived parameters.
# a: Namespace, output of GetArgs()
def GetPar(a):
    d = dict(vars(a))

    we = a.we
    la = a.la
    nx = a.nx
    np = a.np
    brh = a.brh
    vel0 = a.vel0 if a.dim == 3 else [a.vel0[0], a.vel0[1], 0.]
    dom = a.dom
    #mu = max(a.mu0, a.mu0 * a.mu)
    rho = max(a.rho0, a.rho0 * a.rho)
    tmax = a.tmax

    # mesh step
    hx = dom / nx
    # bubble radius
    brx = brh * hx
    bry = brx * a.bryk
    brz = brx
    br = [brx, bry, brz]
    # bubble diameter
    bd = brx * 2

    if a.mode == "vel":
        # velocity magnitude brx * 2
        vm = brx * 2
        sig = vm ** 2 * bd * rho / we
    elif a.mode == "osc":
        # oscillation frequency: sqrt(7) / pi ~= 0.84
        sig = rho * bd ** 3

    # viscosity
    mu = (sig * rho * bd / la) ** 0.5
    # velocity magnitude
    vm = (we * sig / (bd * rho)) ** 0.5
    # velocity
    vel = [x * vm / norm(vel0) for x in vel0]

    assert IsClose(we, rho * norm(vel) ** 2 * bd / sig)
    assert IsClose(la, rho * sig * bd / mu ** 2)

    if a.rho > 1.: # droplet, phase 2 heavier
        mu0 = mu / a.mu
    else: # bubble
        mu0 = mu

    # dump interval
    dumpdt = tmax / a.nfr

    r3 = range(3)
    # bubble center
    bc = [r + hx * o + oo * r for r,o,oo in zip(br,a.bcoh,a.bcor)]
    # second bubble
    b2r = [br[i] * a.b2rr for i in r3]
    b2c = [bc[i] + br[i] * a.b2xr[i] for i in r3]

    if a.pos == "center":
        dc = [a.dom * 0.5 for i in r3] # domain center
        c = [(bc[i] + b2c[i]) * 0.5 for i in r3]
        bc = [bc[i] - c[i] + dc[i] for i in r3]
        b2c = [b2c[i] - c[i] + dc[i] for i in r3]

    # derived variables to return
    kk = ['sig', 'vel', 'dumpdt', 'br', 'bc', 'b2r', 'b2c', 'mu0']
    for k in kk:
        d[k] = locals()[k]

    return d

# Returns text for b.dat
# cc: center, shape (n,3)
# rr: size, shape (n,3)
def GetBubText(cc, rr):
    t = ""
    for c,r in zip(cc,rr):
        t +=  "{:} {:} {:} {:} {:} {:}\n".format(*(c + r))
    return t

def GetBubName():
    return 'b.dat'

def WriteBub(cc, rr):
    t = GetBubText(cc, rr)
    n = GetBubName()
    open(n, 'w').write(t)

def GetParName():
    return "par.py"

def WritePar(t):
    n = GetParName()
    open(n, 'w').write(t)

# Returns python code defining config derived from c.
# c: Namespace, output of GetArgs()
def GetParPy(c):
    return GetDictPy(GetConf(c))

# Writes file with dimension
def WriteDim(dim):
    open("dim", 'w').write(str(dim))
