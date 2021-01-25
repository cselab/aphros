import numpy as np
import os
import re
import scipy
import scipy.interpolate

def read_vtk(fn):
    s = 0    # initial state
    n = None # num points
    nc = None # num cells
    vn = None # current variable name
    d = dict()
    i = None # current array index
    x = None; y = None; z = None;  # point arrays
    a = None # current array
    cc = None  # cells
    ct = None  # cell types
    m = None # temporary buffer
    with open(fn) as f:
        for l in f:
            if s == 0: # parse header
                if l.count("POINTS"):
                    n = int(re.findall("\D*(\d*)\D*", l)[0])
                    assert n > 0
                    i = 0
                    x = np.empty(n); y = np.empty(n); z = np.empty(n)
                    s = 10
                elif l.count("POINT_DATA"):
                    nd = int(re.findall("\D*(\d*)", l)[0])
                    assert nd == n
                    s = 20
                elif l.count("CELLS"):
                    sp = l.split()
                    m = int(sp[1])
                    nc = m if nc is None else nc
                    assert nc == m
                    assert nc > 0
                    sz = int(sp[2])
                    assert (sz - nc) % nc == 0
                    shape = (nc, (sz - nc) // nc) # assume same size
                    cc = np.empty(shape, dtype=int)
                    i = 0
                    s = 30
                elif l.count("CELL_TYPES"):
                    m = int(l.split()[1])
                    nc = m if nc is None else nc
                    assert nc == m
                    assert nc > 0
                    ct = np.empty(nc, dtype=int)
                    i = 0
                    s = 40
            elif s == 10: # read point coordinates
                sp = l.split()
                if len(sp) == 0:
                    continue
                assert len(sp) == 3
                x[i], y[i], z[i] = map(float, sp)
                i += 1
                if i >= n:
                    d['x'] = x; d['y'] = y; d['z'] = z
                    s = 0
            elif s == 20: # read points header
                sp = l.split()
                if len(sp) == 0:
                    continue
                assert sp[0] == "SCALARS"
                assert sp[2] == "float"
                vn = sp[1]
                a = np.empty(n)
                i = 0
                s = 21
            elif s == 21: # skip next line
                sp = l.split()
                assert sp[0] == "LOOKUP_TABLE"
                assert sp[1] == "default"
                s = 22
            elif s == 22: # read point data
                a[i] = float(l)
                i += 1
                if i >= n:
                    d[vn] = a
                    s = 20
            elif s == 30: # read cell point lists
                sp = l.split()
                if len(sp) == 0:
                    continue
                assert len(sp) == int(sp[0]) + 1
                cc[i] = np.array(sp[1:], dtype=int)
                i += 1
                if i >= nc:
                    d['cells'] = cc
                    s = 0
            elif s == 40: # read cell types
                sp = l.split()
                if len(sp) == 0:
                    continue
                assert len(sp) == 1
                ct[i] = int(sp[0])
                i += 1
                if i >= nc:
                    d['cell_types'] = ct
                    s = 0
    return d

# Interpolate points to uniform grid assuming square mesh
# xp, yp -- coords
# fields -- list of data arrays
# n=sqrt(len(xp)) -- mesh size
# Return: (x1, y1, r)
# x1, y1 -- coords 1D
# r -- list of data arrays (order: ...) TODO: order
def interp_to_uniform(xp, yp, fields=[], n=None):
    if n is None:
        n = int(len(xp) ** 0.5 + 0.5)

    x1 = np.linspace(xp.min(), xp.max(), n)
    y1 = np.linspace(yp.min(), yp.max(), n)

    x, y = np.meshgrid(x1, y1)
    pts = np.vstack((xp, yp)).T

    r = []
    for fp in fields:
        f = scipy.interpolate.griddata(pts, fp, (x, y), method="nearest")
        assert np.all(np.isfinite(f))
        r.append(f)

    return x1, y1, r

