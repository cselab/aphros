import os
import re
from aphros.vtk import ReadVtkPoly
try:
    import numpy as np
except ImportError:
    pass


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
        s = np.array(ll[0].split(), dtype=int)
        # shape z,y,x
        ss = tuple(reversed(s))
        # data flat
        u = np.array(ll[1].split(), dtype=float)
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
        '<DataItem.*Dimensions="(\d*) (\d*) (\d*)".*?> *([a-z0-9_.]*)',
        text)[0]
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


def write_raw_xmf(xmfpath,
                  rawpath,
                  count,
                  spacing=(1, 1, 1),
                  name='data',
                  precision=8):
    '''
    Writes XMF metadata for a `.raw` datafile.
    xmfpath: path to output `.xmf` file
    rawpath: path to binary `.raw` file to be linked
    count: array size as (Nz, Ny, Nx)
    name: name of field
    '''

    txt = '''\
<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
   <Grid Name="mesh" GridType="Uniform">
     <Topology TopologyType="{dim}DCORECTMesh" Dimensions="{nodes*}"/>
     <Geometry GeometryType="{geomtype}">
       <DataItem Name="Origin" Dimensions="{dim}" NumberType="Float" Precision="8" Format="XML">
         {origin*}
       </DataItem>
       <DataItem Name="Spacing" Dimensions="{dim}" NumberType="Float" Precision="8" Format="XML">
         {spacing*}
       </DataItem>
     </Geometry>
     <Attribute Name="{name}" AttributeType="Scalar" Center="Cell">
       <DataItem ItemType="HyperSlab" Dimensions="{countd*}" Type="HyperSlab">
           <DataItem Dimensions="3 {dim}" Format="XML">
             {start*}
             {stride*}
             {count*}
           </DataItem>
           <DataItem Dimensions="{bindim*}" Seek="{seek}" Precision="{precision}" NumberType="{type}" Format="Binary">
             {binpath}
           </DataItem>
       </DataItem>
     </Attribute>
   </Grid>
 </Domain>
</Xdmf>
'''

    def tostrrev(v):
        return ' '.join(map(str, reversed(v)))
    def tostr(v):
        return ' '.join(map(str, v))

    info = dict()
    dim = 3
    info['name'] = name
    info['dim'] = dim
    info['origin'] = tostrrev([0] * dim)
    info['spacing'] = tostrrev(spacing)
    info['start'] = tostrrev([0] * dim)
    info['stride'] = tostrrev([1] * dim)
    info['count'] = tostr(count)
    info['bindim'] = tostr(count)
    info['countd'] = tostr(count)
    info['nodes'] = tostr([a + 1 for a in count])
    info['precision'] = precision
    if precision == 8:
        info['type'] = 'Double'
    else:
        info['type'] = 'Float'
    info['binpath'] = rawpath
    info['seek'] = '0'
    info['geomtype'] = 'ORIGIN_DXDYDZ'
    # Remove '*' which are only used in `aphros/src/dump/xmf.ipp`.
    txt = txt.replace('*}', '}')
    txt = txt.format(**info)

    with open(xmfpath, 'w') as f:
        f.write(txt)


def write_raw_with_xmf(u,
                       xmfpath,
                       rawpath=None,
                       spacing=(1, 1, 1),
                       name='data'):
    '''
    Writes binary data in raw format with XMF metadata.
    u: np.ndarray to write, shape (Nz, Ny, Nx)
    spacing: cell size
    name: name of field
    '''
    if len(u.shape) != 3:
        u = u.reshape((1, ) + u.shape)
    if len(spacing) != 3:
        spacing = list(spacing) + [min(spacing)]
    precision = 8 if u.dtype == np.float64 else 4
    if rawpath is None:
        rawpath = os.path.splitext(xmfpath)[0] + ".raw"
    write_raw_xmf(xmfpath, rawpath, u.shape, spacing, name, precision)
    u.tofile(rawpath)
    return xmfpath


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
