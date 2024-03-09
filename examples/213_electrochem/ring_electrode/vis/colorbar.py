#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse


g_colormap_names = [
    "rainbow",
    "coolwarm",
    "yellow",
    "geo",
]

def get_colormap_data(name):
    data = None
    assert name in g_colormap_names
    if name == "rainbow":
        data = [
            0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.143, 0.0,
            0.0, 0.360784313725, 0.285, 0.0, 1.0, 1.0, 0.429, 0.0,
            0.501960784314, 0.0, 0.571, 1.0, 1.0, 0.0, 0.714, 1.0,
            0.380392156863, 0.0, 0.857, 0.419607843137, 0.0, 0.0, 1.0,
            0.878431372549, 0.301960784314, 0.301960784314
        ]
    elif name == "coolwarm":
        data = [
            0.0, 0.0, 0.0, 0.34902, 0.03125000000000003, 0.039216, 0.062745,
            0.380392, 0.06250000000000006, 0.062745, 0.117647, 0.411765,
            0.09374999999999994, 0.090196, 0.184314, 0.45098,
            0.12499999999999997, 0.12549, 0.262745, 0.501961, 0.15625,
            0.160784, 0.337255, 0.541176, 0.18750000000000003, 0.2, 0.396078,
            0.568627, 0.21875000000000006, 0.239216, 0.454902, 0.6,
            0.24999999999999994, 0.286275, 0.521569, 0.65098,
            0.28124999999999994, 0.337255, 0.592157, 0.701961, 0.3125,
            0.388235, 0.654902, 0.74902, 0.34375, 0.466667, 0.737255, 0.819608,
            0.37500000000000006, 0.572549, 0.819608, 0.878431, 0.40625,
            0.654902, 0.866667, 0.909804, 0.4375, 0.752941, 0.917647, 0.941176,
            0.46875, 0.823529, 0.956863, 0.968627, 0.5, 0.988235, 0.960784,
            0.901961, 0.5, 0.941176, 0.984314, 0.988235, 0.52, 0.988235,
            0.945098, 0.85098, 0.54, 0.980392, 0.898039, 0.784314, 0.5625,
            0.968627, 0.835294, 0.698039, 0.59375, 0.94902, 0.733333, 0.588235,
            0.625, 0.929412, 0.65098, 0.509804, 0.65625, 0.909804, 0.564706,
            0.435294, 0.6875, 0.878431, 0.458824, 0.352941, 0.71875, 0.839216,
            0.388235, 0.286275, 0.7500000000000001, 0.760784, 0.294118,
            0.211765, 0.78125, 0.701961, 0.211765, 0.168627, 0.8125, 0.65098,
            0.156863, 0.129412, 0.84375, 0.6, 0.094118, 0.094118, 0.875,
            0.54902, 0.066667, 0.098039, 0.9062500000000001, 0.501961, 0.05098,
            0.12549, 0.9375, 0.45098, 0.054902, 0.172549, 0.96875, 0.4,
            0.054902, 0.192157, 1.0, 0.34902, 0.070588, 0.211765
        ]
    elif name == "geo":
        data = [
            0.0, 0.0, 0.6039215686274509, 0.8705882352941177, 0.5, 1.0, 1.0,
            1.0, 1.0, 1.0, 0.12156862745098039, 0.3568627450980392
        ]
    elif name == "yellow":
        data = [
            0.0, 1.0, 1.0, 0.988235, 0.002, 1.0, 1.0, 0.988235,
            0.05000000000000001, 0.984314, 0.988235, 0.843137,
            0.10000000000000002, 0.988235, 0.988235, 0.741176, 0.15, 0.980392,
            0.968627, 0.654902, 0.20000000000000004, 0.980392, 0.945098,
            0.576471, 0.25, 0.968627, 0.905882, 0.486275, 0.3, 0.968627,
            0.862745, 0.388235, 0.3499999999999999, 0.960784, 0.803922,
            0.286275, 0.4000000000000001, 0.94902, 0.741176, 0.219608, 0.45,
            0.941176, 0.678431, 0.14902, 0.5, 0.929412, 0.607843, 0.094118,
            0.55, 0.921569, 0.545098, 0.054902, 0.6, 0.909804, 0.486275,
            0.035294, 0.65, 0.890196, 0.411765, 0.019608, 0.6999999999999998,
            0.8, 0.305882, 0.0, 0.7500000000000001, 0.760784, 0.239216, 0.0,
            0.8000000000000002, 0.678431, 0.180392, 0.011765, 0.85, 0.6,
            0.121569, 0.023529, 0.9, 0.501961, 0.054902, 0.031373, 0.95, 0.4,
            0.039216, 0.058824, 1.0, 0.301961, 0.047059, 0.090196
        ]
    else:
        assert False, "Unknown colormap=" + name
    data = np.reshape(data, (-1, 4))
    return data


def get_cmap(name):
    data = get_colormap_data(name)
    nodes = data[:, 0]
    colors = data[:, 1:]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        name, list(zip(nodes, colors)))
    return cmap


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('name',
                        type=str,
                        choices=g_colormap_names,
                        default="rainbow",
                        help="Name of colormap")
    parser.add_argument('--output', type=str, default="cbar.pdf")
    parser.add_argument('--size', nargs=2, type=float, default=[3, 1])
    parser.add_argument('--vertical', action="store_true")
    parser.add_argument('--horizontal', action="store_true")
    parser.add_argument('--dpi', type=float, default=300)
    parser.add_argument('--vmin', type=float, default=0)
    parser.add_argument('--vmax', type=float, default=1)
    parser.add_argument('--ticks', type=int, default=5)
    parser.add_argument('--linewidth', type=float, default=0.5)
    parser.add_argument('--levels', type=int, default=20)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': 'cmr10',
        'mathtext.fontset': 'cm',
        'axes.unicode_minus': False,
        'axes.linewidth': args.linewidth,
        'axes.linewidth': args.linewidth,
        'mathtext.rm': 'serif',
    })

    cmap = get_cmap(args.name)

    vertical = not args.horizontal

    if vertical:
        args.size = list(sorted(args.size))
    else: # Horizontal.
        args.size = list(reversed(sorted(args.size)))

    plt.figure(figsize=args.size)

    vmin = args.vmin
    vmax = args.vmax
    if args.levels:
        plt.contourf([[vmin, vmax]] * 2,
                     levels=args.levels - 1,
                     vmin=vmin,
                     vmax=vmax,
                     cmap=cmap)
    else: # Continuous.
        plt.imshow([[vmin, vmax]], cmap=cmap)
    plt.gca().set_visible(False)


    cax = plt.axes([0, 0, 0.3, 1] if vertical else [0, 0, 1, 0.2])
    ticks = np.linspace(vmin, vmax, args.ticks)
    cbar = plt.colorbar(orientation="vertical" if vertical else "horizontal",
                        cax=cax,
                        ticks=ticks)
    cbar.ax.tick_params(width = args.linewidth)
    print(args.output)
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')


if __name__ == "__main__":
    main()
