#!/usr/bin/env python

import numpy as np
import re
import sys
import inspect
import os

def WE(m):
    sys.stderr.write(str(m) + "\n")

def ReadVtkPoly(fn):
    def Assert(cond):
        if not cond:
            caller = inspect.getframeinfo(inspect.stack()[1][0])
            lines = "\n".join(caller[3]).strip()
            filename = os.path.basename(caller.filename)
            lineno = caller.lineno
            WE("\n{:}:{:} {:}".format(filename, lineno, lines))
            # FIXME: not line number but number of iterations
            WE("Failing at input file line {:} in state s={:}:\n{:}".format(
                lnum + 1, s, l.strip()))
            exit(1)
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

    points = None
    poly = None
    dim = 3
    num_points = None
    num_poly = None
    cell_fields = []
    cell_field_name = None
    with open(fn) as f:
        for lnum,l in enumerate(f):
            if not l.strip():
                continue
            if s == 0: # check header
                Assert("# vtk" in l)
                s += 1
            elif s == 1: # skip comment
                s += 1
            elif s == 2:
                Assert("ASCII" in l)
                s += 1
            elif s == 3:
                Assert("DATASET POLYDATA" in l)
                s += 1
            elif s == 4:
                Assert("POINTS" in l)
                num_points = int(re.findall("\D*(\d*)\D*", l)[0])
                points = np.empty((num_points, dim))
                points = np.loadtxt(f, max_rows=num_points)
                Assert(points.shape[0] == num_points)
                WE("Read {:} points".format(points.shape[0]))
                s += 1
            elif s == 5:
                Assert("POLYGONS" in l)
                m = re.findall("\D*(\d*)\s*(\d*)", l)[0]
                num_poly = int(m[0])
                poly = np.loadtxt(f, max_rows=num_poly)
                Assert(poly.shape[0] == num_poly)
                WE("Read {:} polygons".format(poly.shape[0]))
                s += 1
            elif s == 6:
                Assert("CELL_DATA" in l)
                n = int(re.findall("\D*(\d*)", l)[0])
                Assert(n == num_poly)
                s = 20
            elif s == 20: # read cell field
                Assert("SCALARS" in l)
                cell_field_name = re.findall("SCALARS\s*(\S+)", l)[0]
                s += 1
            elif s == 21:
                Assert("LOOKUP_TABLE" in l)
                s += 1
            elif s == 22:
                u = np.loadtxt(f, max_rows=num_poly-1)
                Assert(u.shape[0] == num_poly)
                print("Read cell field '{:}'".format(cell_field_name))
                cell_fields.append(u)
                s += 1

                if "POINTS" in l:
                    n = int(re.findall("\D*(\d*)\D*", l)[0])
                    assert n > 0
                    i = 0
                    x = np.empty(n); y = np.empty(n); z = np.empty(n)
                    s = 10
                elif "POINT_DATA" in l:
                    nd = int(re.findall("\D*(\d*)", l)[0])
                    assert nd == n
                    s = 20
                elif "CELLS" in l:
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
                elif "CELL_TYPES" in l:
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
