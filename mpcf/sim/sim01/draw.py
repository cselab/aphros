#!/usr/bin/env python

import h5py
import glob
import os
import numpy as np
from matplotlib import use as matplotlib_use
matplotlib_use("Agg")
import matplotlib.pyplot as plt
import re

# Read array from given data file
# axis=0: y
# axis=1: x
# e.g. u[y,x]
def GetArray(data_path, component=None, coarse_size=None, radius=1., dataset="data"):
    f = h5py.File(data_path, 'r')
    dset = f[dataset]
    u = np.array(dset[()])
    f.close()
    if radius != 1.:
        a = 0.5 * (1. - radius)
        b = 0.5 * (1. + radius)
        s = u.shape
        u = u[int(a * s[0]):int(b * s[0]), int(a * s[1]):int(b * s[1])]
    if coarse_size is not None:
        skip = [max(1, s / coarse_size) for s in u.shape]
        u = u[::skip[0], ::skip[1], :]

    if u.shape[2] == 1:
        return u[:,:,0]
    elif component is not None:
        return u[:,:,component]
    else:
        return u

def GetGradient(u, dx=1.):
    gx = 0.5 * (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / dx
    gy = 0.5 * (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / dx
    return gx, gy

def GetMagnitude(gx, gy):
    return np.sqrt(np.square(gx) + np.square(gy))

def GetSchlierenFromMagnitude(magnitude, upper=None, k=1.):
    #upper = np.percentile(magnitude, 99)
    if upper is None:
        upper = np.max(magnitude)
    print("upper=", upper)
    return np.exp(-k * magnitude / upper)

# Make contour and line plots
def DrawArray(u_arg, basename, vmin=None, vmax=None, text="", schlupper=None,
        contour_scale=1., line_scale=1., data_alpha=None):
    u = np.array(u_arg)
    if vmin is None:
        vmin = u.min()
    if vmax is None:
        vmax = u.max()

    ptp = vmax - vmin

    # create contour figure
    fig = plt.figure(figsize=(5.,5.))
    fig.patch.set_facecolor('black')
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    im = ax.imshow(np.flipud(u), vmin=vmin, vmax=vmax, cmap="coolwarm")
    if data_alpha is not None:
        ax.contour(np.flipud(data_alpha), [0.5], colors='k', linewidths=0.1)
    # colorbar
    if 0:
        cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.7])
        ticks = [vmin, 0.5 * (vmin + vmax), vmax]
        ticklabels = [vmin, 0.5 * (vmin + vmax), vmax]
        cbar = fig.colorbar(im, ax=ax, cax=cbaxes,
                            ticks=ticks, orientation='vertical')
        cbar.ax.set_yticklabels(ticklabels, color='w')
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
        plt.text(0.97, 0.97,text,
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform = fig.transFigure, color='w')
    contour_path = "contour_{:}.png".format(basename)
    print("Save to {:}".format(contour_path))
    plt.savefig(contour_path, dpi=1080/fig.get_size_inches()[0], facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()

    # create schlieren figure
    fig = plt.figure(figsize=(5.,5.))
    fig.patch.set_facecolor('black')
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    gx, gy = GetGradient(u)
    mag = GetMagnitude(gx, gy)
    us = GetSchlierenFromMagnitude(mag, k=10., upper=schlupper)
    im = ax.imshow(np.flipud(us), cmap='gray', vmin=0., vmax=1.)
    if data_alpha is not None:
        ax.contour(np.flipud(data_alpha), [0.5], colors='k', linewidths=0.1)
    plt.text(0.97, 0.97,text,
             horizontalalignment='right',
             verticalalignment='top',
             transform = fig.transFigure, color='w')
    schl_path = "schl_{:}.png".format(basename)
    print("Save to {:}".format(schl_path))
    plt.savefig(schl_path, dpi=1080/fig.get_size_inches()[0], facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()

    # create line figure
    fig = plt.figure()
    uc = u[u.shape[0] // 2]
    x = np.linspace(-1., 1., u.shape[1])
    plt.plot(x, uc)
    linepad = 0.5
    plt.ylim(vmin - linepad * ptp, vmax + linepad * ptp)
    plt.grid()
    plt.tight_layout()
    line_path = "line_{:}.png".format(basename)
    print("Save to {:}".format(line_path))
    plt.savefig(line_path)
    plt.close()


sim = np.load("../sim.npy").item()
ampfac = sim["amp__pinf"]
amb = sim["pinf"]
print("ampfac={:}, amb={:}".format(ampfac, amb))


files = sorted(glob.glob("data*.h5"))[::1]
for idx,f in enumerate(files):
  name = os.path.basename(f)
  print("{:}/{:} {:}".format(idx + 1, len(files), f))
  u = GetArray(f)
  basename  = os.path.splitext(os.path.basename(f))[0]
  if len(re.findall("-p_", name)):
    alpha_path = "../alpha/" + f.replace("-p_", "-a2_")
    data_alpha = GetArray(alpha_path)
    schlupper = 60000. * ampfac
    DrawArray(u, basename, vmin=amb * (1. - ampfac), vmax=amb * (1. + ampfac),
            schlupper=schlupper, data_alpha=data_alpha)
  elif len(re.findall("-a2_", name)):
    DrawArray(u, basename, vmin=0., vmax=1.)
  else:
    print("Warning: field name not recognized. Use auto range")
    DrawArray(u, basename)

