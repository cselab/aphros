#!/usr/bin/env python3

import numpy as np
import re
import sys
import inspect
import os

def WE(m):
    sys.stderr.write(str(m) + "\n")

def ReadVtkPoly(fn, verb=False):
    def Assert(cond, msg=""):
        if not cond:
            caller = inspect.getframeinfo(inspect.stack()[1][0])
            lines = "\n".join(caller[3]).strip()
            filename = os.path.basename(caller.filename)
            lineno = caller.lineno
            WE("\n{:}:{:} {:}".format(filename, lineno, lines))
            WE("Failing at iteration {:} in state s={:}".format(lnum + 1, s))
            if msg: WE(str(msg))
            WE("Current input line:\n{:}".format(l.strip()))
            WE("Next line would be:\n{:}".format(f.readline().strip()))
            exit(1)
    s = 0    # initial state

    points = None
    poly = None
    dim = 3
    num_points = None
    num_poly = None
    cell_fields = dict()
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
                Assert(points.shape[1] == 3)
                if verb: WE("Read {:} points".format(points.shape[0]))
                s += 1
            elif s == 5:
                Assert("POLYGONS" in l)
                m = re.findall("\D*(\d*)\s*(\d*)", l)[0]
                num_poly = int(m[0])
                poly = np.loadtxt(f, max_rows=num_poly)
                Assert(poly.shape[0] == num_poly)
                if verb: WE("Read {:} polygons".format(poly.shape[0]))
                s += 1
            elif s == 6:
                if "CELL_DATA" in l:
                    n = int(re.findall("\D*(\d*)", l)[0])
                    Assert(n == num_poly)
                    s = 20
                elif "POINT_DATA" in l:
                    pass
            elif s == 20: # read cell field
                if "SCALARS" in l:
                    cell_field_name = re.findall("SCALARS\s*(\S+)", l)[0]
                    s += 1
                else:
                    s = 6
            elif s == 21:
                Assert("LOOKUP_TABLE" in l)
                u = np.loadtxt(f, max_rows=num_poly)
                Assert(u.shape[0] == num_poly, ["u.shape=", u.shape])
                if verb: WE("Read cell field '{:}'".format(cell_field_name))
                cell_fields[cell_field_name] = u
                s = 20
    return points, poly, cell_fields
