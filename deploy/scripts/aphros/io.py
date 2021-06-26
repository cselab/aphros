import numpy
import os
import re
import numpy as np

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
