#!/usr/bin/env python3

import numpy as np
import re
import sys
import argparse

import aphros
import plottools
import matplotlib.pyplot as plt
aphros.plot.ApplyParams(plt)

def printerr(m):
    sys.stderr.write('{:}\n'.format(m))
    sys.stderr.flush()


def get_polar(u):
    '''
    Returns profile versus angle relative to the domain center ignoring NaN.
    u: 2D array, shape (ny, nx)
    '''
    iy, ix = np.where(np.isfinite(u))
    shape = u.shape
    cx = (shape[0] - 1) * 0.5
    cy = (shape[1] - 1) * 0.5
    dx = np.array(ix) - cx
    dy = np.array(iy) - cy
    angle = np.arctan2(dx, dy)
    argsort = np.argsort(angle)
    return angle[argsort], u[iy, ix][argsort]

def parse_xmf(xmfpath):
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
    return shape, rawpath

def read_raw(xmfpath):
    '''
    Returns array from scalar field in raw format.
    xmfpath: path to xmf metadata file
    '''
    shape, rawpath = parse_xmf(xmfpath)
    u = np.fromfile(rawpath).reshape(shape)
    return u

def parse_args():
    parser = argparse.ArgumentParser(description="test")
    parser.add_argument('--xmfpath',
                        type=str,
                        default="k_0000.xmf",
                        help="Path to xmf file with curvature field")
    parser.add_argument('--output',
                        type=str,
                        default="curvature.pdf",
                        help="Path to output image")
    parser.add_argument('--radius',
                        type=float,
                        default="1",
                        help="Radius of circle to compute exact curvature")
    return parser.parse_args()



if __name__ == "__main__":
    args = parse_args()
    k = read_raw(args.xmfpath)
    angle, k = get_polar(k[0, :, :])

    k_exact = 1 / args.radius
    k_error = (k - k_exact) / k_exact

    fig, ax = plt.subplots()
    ax.plot(angle, k_error, lw=1)
    ax.set_xlabel('angle')
    ax.set_ylabel('curvature error')
    plottools.adjust_ticks(ax)
    printerr(args.output)
    fig.tight_layout()
    fig.savefig(args.output, dpi=300)
