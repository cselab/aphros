#!/usr/bin/env python3

import argparse
import math
import os
import sys


# Returns dictionary with standard parameters {k : (v,h)}
# k: name
# v: default value or type [int, float, list] (then default is None)
# h: description for help
def StdParAll():
    Re = [float, "Reynolds number"]
    Red = [float, "Reynolds number by bubble diameter"]
    Reh = [float, "Reynolds number by channel height"]
    Ca = [float, "capillary number"]
    We = [float, "Weber number"]
    La = [float, "Laplace number"]
    Ga = [float, "Galilei number"]
    Fr = [float, "Froude number (zero for no gravity)"]
    Ga = [float, "Galilei number"]
    Eo = [float, "Eotvos number"]
    Oh = [float, "Ohnesorge number"]
    np = [1, "number of processors"]
    nx = [int, "mesh size"]
    cpr = [float, "cells per bubble radius"]
    br = [0.125, "bubble radius relative to dom"]
    tmax = [1., "total time"]
    nfr = [100, "number of frames"]
    mu0 = [1., "dynamic viscosity of carrier fluid"]
    mu = [1., "dynamic viscosity of bubble relative to carrier"]
    rho0 = [1., "density of carrier fluid"]
    rho = [1., "density of bubble relative to carrier"]
    dim = [3, "dimension, 2 or 3"]
    chsm = [1, "ch smooth steps [rhor>10: 1, rhor>100: 2]"]
    gesm = [1, "ge smooth steps"]
    wallx = [0, "wall in x direction, else periodic"]
    wally = [0, "wall in y direction, else periodic"]
    wallz = [0, "wall in z direction, else periodic"]
    vel0 = [[1., 0.8, 0.6], "direction of velocity"]
    g0 = [[0., -1., 0.], "direction of gravity"]
    bcoh = [[0., 0., 0.], "offset of bubble center relative to h"]
    bcod = [[0., 0., 0.], "offset of bubble center relative to dom"]
    bcor = [[0., 0., 0.], "offset of bubble center relative to br"]
    bryk = [1., "stretching factor for bubble size in y"]
    dom = [1., "domain size"]
    ch = [1, "run ch", [0, 1]]
    ge = [1, "run ge", [0, 1]]
    pos = ["corner", "initial position of bubbles", ["center", "corner"]]
    mode = [
        "vel", '''timescale; vel: velocity magnitude D;
osc: oscillation frequency sqrt[7] / pi;
mu: viscosity mu = rho * D ^ 2, pois: velocity magnitude dom''',
        ['vel', 'osc', 'mu', 'pois']
    ]
    pois = [0, "Poiseuille profile by force along x", [0, 1]]
    b2rr = [0., "second bubble radius relative to br"]
    b2xr = [[2., 0., 0.], "second bubble displacement realative to br"]
    out = [str, "directory for task and output"]
    meshvel = [0, "moving mesh in ch", [0, 1]]
    nb = [1, "number of bubbles"]
    symm = [0, "symmetry conditions at y=0.5 and z=0.5"]
    symmxy = [0, "symmetry conditions at x=0.5 and y=0.5"]
    symmcorn = [0, "symmetry conditions at x,y,z=0 corner"]
    dumpdefault = [0, "add dumpformat=default to ch"]
    part_np = [7, "particles in one string"]
    part_ns = [3, "strings per cell"]
    part_hp = [4., "string length"]
    part_circ = [0., "segcirc factor"]
    part_eps = [1e-5, "convergence tolerance"]
    part_eta = [0.5, "relaxation factor"]
    part_itermax = [20, "maximum number of iterations"]
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
    return rho * vel**2 * x / sig


def GetLa(x, rho, mu, sig):
    return rho * sig * x / mu**2


def GetGa(x, rho, mu, sig):
    return rho * sig * x / mu**2


# Returns a selection from StdParAll().
# kk: list of names
def GetStdPar(*kk):
    if type(kk) == str:
        kk = [kk]
    d = StdParAll()
    r = dict()  # result
    for k in kk:
        r[k] = d[k]
    return r


def IsClose(a, b):
    return abs(a - b) < 1e-10


# Returns Namespace with parameters from c overriden by args
# cc: dict with parameters as in StdParAll()
# desc: description for command line help
def GetArgs(cc, desc=None):
    p = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    for k in sorted(cc):
        c = cc[k]
        v = c[0]  # default value
        h = c[1]  # help message
        o = c[2] if len(c) > 2 else None  # choices
        if v in [int, float, list, str]:
            t = v
            v = None
        else:
            t = type(v)
            v = v

        n = '-' + k

        if t == list:
            m = ('X', 'Y', 'Z')
            assert v is None or len(v) == len(m)
            t = type(v[0]) if v is not None else float
            p.add_argument(n,
                           nargs=len(m),
                           choices=o,
                           metavar=m,
                           default=v,
                           help=h,
                           type=t)
        else:
            p.add_argument(n, nargs='?', choices=o, default=v, help=h, type=t)

    a = p.parse_args()
    return a


def norm(v):
    assert len(v) == 3
    return sum([a**2 for a in v])**0.5


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
        sig = vm**2 * bd * rho / we
    elif a.mode == "osc":
        # oscillation frequency: sqrt(7) / pi ~= 0.84
        sig = rho * bd**3

    # viscosity
    mu = (sig * rho * bd / la)**0.5
    # velocity magnitude
    vm = (we * sig / (bd * rho))**0.5
    # velocity
    vel = [x * vm / norm(vel0) for x in vel0]

    assert IsClose(we, rho * norm(vel)**2 * bd / sig)
    assert IsClose(la, rho * sig * bd / mu**2)

    if a.rho > 1.:  # droplet, phase 2 heavier
        mu0 = mu / a.mu
    else:  # bubble
        mu0 = mu

    # dump interval
    dumpdt = tmax / a.nfr

    r3 = range(3)
    # bubble center
    bc = [r + hx * o + oo * r for r, o, oo in zip(br, a.bcoh, a.bcor)]
    # second bubble
    b2r = [br[i] * a.b2rr for i in r3]
    b2c = [bc[i] + br[i] * a.b2xr[i] for i in r3]

    if a.pos == "center":
        dc = [a.dom * 0.5 for i in r3]  # domain center
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
    for c, r in zip(cc, rr):
        t += "{:} {:} {:} {:} {:} {:}\n".format(*(c + r))
    return t


def GetBubName():
    return 'b.dat'


def WriteBub(cc, rr):
    t = GetBubText(cc, rr)
    n = GetBubName()
    open(n, 'w').write(t)


def GetParName():
    return "par.py"


# Returns python code defining variables from dict c
def GetDictPy(c):
    t = ""
    for k in sorted(c):
        v = c[k]
        v = '"{:}"'.format(v) if type(v) == str else str(v)
        t += "{:} = {:}\n".format(k, v)

    return t


# Writes file with parameters to current directory
# c: dict with parameters
def WritePar(c):
    open(GetParName(), 'w').write(GetDictPy(c))


# Returns dict with parameters
# b: directory containing 'par.conf'
def ReadPar(b):
    f = os.path.join(b, GetParName())
    c = {}
    exec(open(f).read(), None, c)
    return c


# Writes file with dimension
def WriteDim(dim):
    open("dim", 'w').write(str(dim))


# Block size
def GetBs(dim):
    return 8 if dim == 3 else 16


# Returns arguments for base.makefile
# c: dict with parameters, output of ReadPar()
# Required parameters: np, dim, nx
def GetMakeArg(c):
    np = c['np']
    loc = True if np == 1 else False
    dim = c['dim']
    d3 = (dim == 3)

    # minimal size in z
    mz = 1 if loc else 2

    nx = c['nx']
    ny = nx
    nz = nx if d3 else mz

    # block size
    bx = GetBs(dim)
    by = bx
    bz = bx if d3 else mz

    return "m='{nx} {ny} {nz}' bs='{bx} {by} {bz}' np={np}".format(**locals())


# Returns lowest power of 2 which is an upper limit
def Upper2(a):
    r = 0
    while a > 2**r:
        r += 1
    return 2**r


def sh(s, fatal=True, silent=True):
    if not silent:
        print(s)
    r = os.system(s)
    assert not fatal or r == 0
