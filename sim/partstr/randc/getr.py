#!/usr/bin/env python

# ri: radius index
# Returns:
# r: radius
# nx: mesh size
# cpr: r*nx
def GetR(ri):
    cpr = 2. ** (ri * 0.25 - 1)
    nxb = 16   # block size in ch
    nx = max(4, int((2. + cpr) * 2. + 0.5) + nxb - 1) // nxb * nxb
    r = cpr / nx
    return r, nx, cpr
