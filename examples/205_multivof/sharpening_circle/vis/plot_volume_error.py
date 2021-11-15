#!/usr/bin/env python3

import numpy as np
import argparse
from argparse import Namespace
import plottools
import matplotlib.pyplot as plt
import aphros
import os

myname = os.path.splitext(os.path.basename(__file__))[0]


def load_stat(path):
    u = np.genfromtxt(path, names=True)
    res = Namespace(**{k: u[k] for k in u.dtype.names})
    return res


parser = argparse.ArgumentParser()
parser.add_argument('dir', nargs='?', type=str, default='.')
parser.add_argument('--output', type=str, default="volume_error.pdf")
parser.add_argument('--ylim', nargs='*', type=float, default=[-1, 1])
parser.add_argument('--yscpow', nargs='*', type=int, default=-15)
args = parser.parse_args()

plottools.apply_params(plt)

stat = load_stat(os.path.join(args.dir, "stat.dat"))

fig, ax = plt.subplots(figsize=(1.4,1.2))
yscpow = args.yscpow
if args.ylim:
    ax.set_ylim(*args.ylim)
ax.set_xlabel('step')
ax.set_ylabel(r'volume error $[10^{{{:}}}]$'.format(yscpow))
ax.plot(stat.step, stat.vol2_diff * 10 ** (-yscpow))
plottools.savefig(fig, args.output, transparent=True)
