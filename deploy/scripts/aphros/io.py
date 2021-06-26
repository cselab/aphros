import numpy
import os
import re
import numpy as np
from aphros.vtk import ReadVtkPoly

def ReadPlain(fn):
    '''
    Reads uniform grid data.
    p: path

    Format:
    <nx> <ny> <nz>
    <u[0,0,0]> <u[0,0,1]> ...

    Returns:
    array of shape (nx, ny, nz) or None if file not found
    '''
    assert os.path.isfile(fn), "No such file: '{:}'".format(fn)
    with open(fn) as f:
        ll = f.readlines()
        # shape x,y,z
        s = numpy.array(ll[0].split(), dtype=int)
        # shape z,y,x
        ss = tuple(reversed(s))
        # data flat
        u = numpy.array(ll[1].split(), dtype=float)
        # data z,y,x
        u = u.reshape(ss)
        return u

def parse_raw_xmf(xmfpath):
    '''
    Returns shape and path to `.raw` file
    xmfpath: path to `.xmf` metadata file
    '''
    with open(xmfpath) as f:
        text = ''.join(f.read().split('\n'))
    m = re.findall(
        '<Xdmf.*<Attribute.*'
        '<DataItem.*<DataItem.*'
        '<DataItem.*Dimensions="(\d*) (\d*) (\d*)".*?> *([a-z0-9_.]*)', text)[0]
    shape = tuple(map(int, m[:3]))
    rawpath = m[3]
    rawpath = os.path.join(os.path.dirname(xmfpath), rawpath)
    return shape, rawpath

def read_raw(xmfpath):
    '''
    Returns array from scalar field in raw format.
    xmfpath: path to xmf metadata file
    '''
    shape, rawpath = parse_raw_xmf(xmfpath)
    u = np.fromfile(rawpath).reshape(shape)
    return u

def read_lines_vtk(path):
    '''
    Returns array of lines from `sm_*.vtk` generated with `spacedim=2`
    [line0_x, line0_y, line1_x, line1_y, ...]
    where line0_x and line1_y are lists of coordinates
    '''
    points, poly, _ = ReadVtkPoly(path)
    # [[[x0a, y0b], [x0a, y0b]], ...]
    lines = points[poly][:, :, :2]
    # [[[x0a, x0b], [y0a, y0b]], ...]
    lines = np.transpose(lines, (0, 2, 1))
    # [[x0a, x0b], [y0a, y0b], ...]
    lines = lines.reshape((-1, 2))
    return lines

def walk_line(edges, i, j, visited, points):
    '''
    Returns line containing all adjacent segments starting from edge
    `(points[i], points[j])`.
    '''
    if (i, j) in visited:
        return []
    tail = [points[i], points[j]]
    visited.add((i, j))
    visited.add((j, i))
    term = False
    while not term:
        term = True
        for jp in edges[j]:
            if (j, jp) not in visited:
                tail.append(points[jp])
                visited.add((j, jp))
                visited.add((jp, j))
                j = jp
                term = False
                break
    head = []
    term = False
    while not term:
        term = True
        for ip in edges[i]:
            if (i, ip) not in visited:
                head.append(points[ip])
                visited.add((i, ip))
                visited.add((ip, i))
                i = ip
                term = False
                break
    return list(reversed(head)) + tail


def join_lines(points, poly):
    # edges[i] is indices of points adjacent to points[i] through a line segment
    edges = [[] for i in range(len(points))]
    for i, j in poly:
        if j not in edges[i]: edges[i].append(j)
        if i not in edges[j]: edges[j].append(i)
    # visited contains (i,j) if `(points[i], points[j])` is visited
    visited = set()
    res = []
    # start from each edge and join all adjacent line segments
    for i in range(len(points)):
        for j in edges[i]:
            line = walk_line(edges, i, j, visited, points)
            if line:
                res.append(line)
    return res


def read_joined_lines_vtk(path):
    '''
    Returns array of lines from `sm_*.vtk` generated with `spacedim=2`
    joining adjacent line segments.
    [line0_x, line0_y, line1_x, line1_y, ...]
    where line0_x and line1_y are lists of coordinates
    '''
    points, poly, _ = ReadVtkPoly(path)
    joined = join_lines(points, poly)
    lines = []
    for line in joined:
        lines.append([p[0] for p in line])
        lines.append([p[1] for p in line])
    return lines


