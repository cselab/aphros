#!/usr/bin/env python3

try:
    import numpy as np
except ImportError:
    pass
import re
import sys
import inspect
import os


def printerr(m):
    sys.stderr.write(str(m) + "\n")

def ReadVtkPoly(f, verbose=False):
    """
    Reads vtk points, polygons and fields from legacy VTK file.
    f: `str` or file-like
        Path to legacy VTK file or file-like object.
    Returns:
    points: `numpy.ndarray`, shape (num_points, 3)
        List of points.
    poly: `list` [`list` [ `int` ]], shape (num_cells, ...)
        Polygons as lists of indices in `points`.
    cell_fields: `dict` [`str`, `numpy.ndarray`] , shape (num_cells,)
        Cell felds indexed by name. Each field has shape (num_cells,).
    """
    def Assert(cond, msg=""):
        if not cond:
            caller = inspect.getframeinfo(inspect.stack()[1][0])
            lines = "\n".join(caller[3]).strip()
            filename = os.path.basename(caller.filename)
            lineno = caller.lineno
            printerr("\n{:}:{:} {:}".format(filename, lineno, lines))
            printerr("Failing at iteration {:} in state s={:}".format(
                lnum + 1, s))
            if msg: printerr(str(msg))
            printerr("Current input line:\n{:}".format(l.strip()))
            printerr("Next line would be:\n{:}".format(f.readline().strip()))
            exit(1)

    class S:
        header, comment, binary, dataset, points, \
        polygons, cell_data, cell_scalars, cell_field = range(9)

    points = None
    poly = None
    dim = 3
    num_points = None
    num_poly = None
    cell_fields = dict()
    cell_field_name = None
    binary = False

    path = None
    if type(f) is str:
        path = f
        f = open(path, 'rb')
    else:
        pass # expect file-like

    s = S.header
    if f:
        for lnum, l in enumerate(f):
            l = str(l)
            if not l.strip():
                continue
            if s == S.header:  # check header
                Assert("# vtk" in l)
                s = S.comment
            elif s == S.comment:  # skip comment
                s = S.binary
            elif s == S.binary:
                Assert("ASCII" in l or "BINARY" in l)
                binary = "BINARY" in l
                s = S.dataset
            elif s == S.dataset:
                Assert("DATASET POLYDATA" in l)
                s = S.points
            elif s == S.points:
                Assert("POINTS" in l)
                dtype = np.float64 if "double" in l else np.float32
                num_points = int(re.findall("\D*(\d*)\D*", l)[0])
                points = np.empty((num_points, dim))
                if binary:
                    dt = np.dtype('>f4')
                    bytes = f.read(3 * num_points * dt.itemsize)
                    points = np.frombuffer(bytes, dtype=dt)
                    points = points.astype(np.float)
                    f.readline()
                else:
                    points = np.fromfile(f,
                                         dtype=np.float,
                                         count=num_points * 3,
                                         sep=' ')
                points = points.reshape((num_points, 3))
                Assert(points.shape[0] == num_points)
                Assert(points.shape[1] == 3)
                if verbose: printerr("Read {:} points".format(points.shape[0]))
                s = S.polygons
            elif s == S.polygons:
                Assert("POLYGONS" in l)
                m = re.findall("\D*(\d*)\s*(\d*)", l)[0]
                num_poly = int(m[0])
                num_ints = int(m[1])
                if binary:
                    dt = np.dtype('>i')
                    bytes = f.read(num_ints * dt.itemsize)
                    ints = np.frombuffer(bytes, dtype=dt)
                    ints = ints.astype(np.int)
                    f.readline()
                else:
                    ints = np.fromfile(f,
                                       dtype=np.int,
                                       count=num_ints,
                                       sep=' ')
                i = 0
                poly = []
                for ip in range(num_poly):
                    n = ints[i]
                    i += 1
                    poly.append(ints[i:i + n])
                    i += n
                Assert(i == num_ints)
                Assert(len(poly) == num_poly)
                if verbose: printerr("Read {:} polygons".format(len(poly)))
                s = S.cell_data
            elif s == S.cell_data:
                if "CELL_DATA" in l:
                    n = int(re.findall("\D*(\d*)", l)[0])
                    Assert(n == num_poly)
                    s = S.cell_scalars
                elif "POINT_DATA" in l:
                    pass
            elif s == S.cell_scalars:  # read cell field
                if "SCALARS" in l:
                    cell_field_name = re.findall("SCALARS\s*(\S+)", l)[0]
                    s = S.cell_field
                else:
                    s = S.cell_data
            elif s == S.cell_field:
                Assert("LOOKUP_TABLE" in l)
                if binary:
                    dt = np.dtype('>f4')
                    bytes = f.read(num_poly * dt.itemsize)
                    u = np.frombuffer(bytes, dtype=dt)
                    u = u.astype(np.float)
                    f.readline()
                else:
                    u = np.fromfile(f, dtype=np.float, count=num_poly, sep=' ')

                Assert(u.shape[0] == num_poly, ["u.shape=", u.shape])
                if verbose:
                    printerr("Read cell field '{:}'".format(cell_field_name))
                cell_fields[cell_field_name] = u
                s = S.cell_scalars

    if path:
        f.close()

    return points, poly, cell_fields

def WriteVtkPoly(f, points, poly, comment=""):
    """
    Writes polygons to ASCII legacy VTK file.
    f: `str` or file-like
        Path to output legacy VTK file or file-like object.
    points: `numpy.ndarray`, shape (num_points, 3)
        List of 3D points.
    poly: `list` [`list` [ `int` ]], shape (num_cells, ...)
        Polygons as lists of indices in `points`.
    """
    path = None
    if type(f) is str:
        path = f
        f = open(path, 'wb')
    else:
        pass # expect file-like

    def write(data):
        if type(data) is str:
            data = data.encode()
        f.write(data)

    write("# vtk DataFile Version 2.0\n")

    write(comment + '\n')

    write("ASCII\n")

    write("DATASET POLYDATA\n")

    num_points = len(points)
    write("POINTS {:} float\n".format(num_points))
    for x in points:
        write("{:} {:} {:}\n".format(*x))

    num_poly = len(poly)
    num_poly_data = len(poly) + sum([len(p) for p in poly])
    write("POLYGONS {:} {:}\n".format(num_poly, num_poly_data))
    for p in poly:
        write(' '.join(map(str, [len(p)] + p)) + '\n')

    if path:
        f.close()

