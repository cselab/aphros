#!/usr/bin/env python3

# ri: radius index
# Returns:
# r: radius
# nxexp: exponent of mesh size, nx=2**nxexp
# cpr: r*nx
def GetR(ri):
    cpr = 2. ** (ri * 0.25 - 1)
    nxexp = max(4, int((2. + cpr) * 2. + 0.5).bit_length())
    nx = 2. ** nxexp
    r = cpr / nx
    return r, nxexp, cpr
