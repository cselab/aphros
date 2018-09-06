#!/usr/bin/env python

import os

def Default():
    we = 1.     # Weber number
    la = 1000.  # Laplace number
    np = 1      # number of processors
    nx = 64     # mesh size
    brh = 4.    # r/h, cells per bubble radius
    bryk = 1.2  # stretching factor for bubble size in y
    bxr = [4., 3., 2.]  # bubble displacement relative to radius
    bcoh = [0.2, 0.4, 0.6] # offset of bubble center relative to cell size
    nfr = 100    # number of frames
    mu = 0.01      # dynamic viscosity
    rho = 1.     # density
    dom = 1.     # domain size

    return locals().copy()

def IsClose(a, b):
    return abs(a - b) < 1e-10

# Updates c from environ
def ParseEnv(c):
    e = os.environ

    # Updates c with c[k]=f(e[k])
    def E(k, f=float):
        if k in e:
            c[k] = f(e[k])

    stov = lambda v: list(map(float, v.split()))

    E("we")
    E("la")
    E("np", int)
    E("nx", int)
    E("brh")
    E("bryk")
    E("bxr", stov)
    E("bcoh", stov)
    E("nfr", int)
    E("mu")
    E("rho")
    E("dom")

# Returns python code defining variables from dict c
# Returns python code defining config based on c.
# c: dict such that c.keys() includes Default().keys()
def GetDictPy(c):
    return '\n'.join(["{:} = {:}".format(k, c[k]) for k in c])

def norm(v):
    assert len(v) == 3
    return sum([a ** 2 for a in v]) ** 0.5

# Returns dictionary with config derived from c.
# c: dict including Default().keys()
def GetConf(c):
    d = dict(c)

    we = c['we']
    la = c['la']
    nx = c['nx']
    np = c['np']
    brh = c['brh']
    bxr = c['bxr']
    mu = c['mu']
    rho = c['rho']
    dom = c['dom']

    # mesh step
    hx = dom / nx
    # bubble radius
    brx = brh * hx
    bry = brx * c['bryk']
    brz = brx
    br = [brx, bry, brz]
    # bubble diameter
    bd = brx * 2
    # surface tension
    sig = la * mu ** 2 / (rho * bd)
    # velocity magnitude
    vm = (we * sig / (bd * rho)) ** 0.5
    # total time
    tmax = norm(bxr) * brx / vm
    # velocity
    vel = [x * vm / norm(bxr) for x in bxr]
    # bubble center
    bc = [hx * (r / hx + 1 + o) for r,o in zip(br,c['bcoh'])]

    assert IsClose(we, rho * norm(vel) ** 2 * bd / sig)
    assert IsClose(la, rho * sig * bd / mu ** 2)

    # dump interval
    dumpdt = tmax / c['nfr']

    # derived variables to return
    kk = ['tmax', 'sig', 'vel', 'dumpdt', 'br', 'bc']
    for k in kk:
        d[k] = locals()[k]

    return d

# Returns text for b.dat
# c: center, shape 3
# r: size, shape 3
def GetBub(c, r):
    return "{:} {:} {:} {:} {:} {:}".format(*c, *r)

# Returns python code defining config derived from c.
# c: dict such that c.keys() includes Default().keys()
def GetConfPy(c):
    return GetDictPy(GetConf(c))

c = Default()
ParseEnv(c)
d = GetConf(c)
dpy = GetConfPy(d)
open("conf.py", 'w').write(dpy)

b = GetBub(d['bc'], d['br'])
open("b.dat", 'w').write(b)

