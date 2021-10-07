import os
import itertools
import functools
from argparse import Namespace, ArgumentParser
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from cycler import cycler
import numpy as np


def cache_to_file(targetbase, update=False, arg0=False):
    """
    Factory for a decorator that caches the result of
    function and stores it to a target file.

    targetbase: base path to cache file
    update: force cache update
    arg0: append cache name by first argument converted to string

    Example: Creates file "_cache_7.pickle" with `int(7)`.

    @cache_to_file("_cache.pickle", arg0=True)
    def F(a):
        return a
    F(7)
    """

    ext = os.path.splitext(targetbase)[1]
    if ext == '.pickle':
        import pickle

        def load(path):
            with open(path, 'rb') as f:
                print("Loading cache '{}'".format(path))
                return pickle.load(f)

        def save(content, path):
            with open(path, 'wb') as f:
                print("Saving cache '{}'".format(path))
                pickle.dump(content, f, protocol=4)
    elif ext == '.json':
        import json

        def load(path):
            with open(path, 'r') as f:
                print("Loading cache '{}'".format(path))
                return json.load(f)

        def save(content, path):
            with open(path, 'w') as f:
                print("Saving cache '{}'".format(path))
                json.dump(content, f)

    elif ext == '.npy':
        import numpy

        def load(path):
            print("Loading cache '{}'".format(path))
            return numpy.load(path, allow_pickle=True).item()

        def save(content, path):
            print("Saving cache '{}'".format(path))
            numpy.save(path, content, allow_pickle=True)
    else:
        raise ValueError(
            "Unrecognized extension '{}', expected .pickle.".format(ext))

    def decorator(func):
        def inner(*args, **kwargs):
            name, ext = os.path.splitext(targetbase)
            if len(args) and arg0:
                name += '_{:}'.format(args[0])
            target = name + ext
            if not update and os.path.isfile(target):
                return load(target)
            result = func(*args, **kwargs)
            d = os.path.dirname(target)
            if d: os.makedirs(d, exist_ok=True)
            save(result, target)
            return result

        return functools.wraps(func)(inner)

    return decorator


# https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit/tree/master/Colours
geodata_qualitative = [
    "#FF1F5B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E", "#F28522",
    "#A0B1BA", "#A6761D", "#E9002D", "#FFAA00", "#00B000"
]

colorscheme = geodata_qualitative[:7]


def get_params():
    params = {
        'font.size': 7,
        'font.family': 'serif',
        'font.serif': 'cmr10',
        'mathtext.fontset': 'cm',
        'mathtext.rm': 'serif',
        'figure.figsize': (2.3, 1.8),
        'figure.autolayout': True,
        'axes.autolimit_mode': 'round_numbers',
        'axes.xmargin': 0,
        'axes.ymargin': 0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.unicode_minus': False,
        'axes.prop_cycle': cycler(color=colorscheme),
        'lines.markersize': 3,
        'lines.linewidth': 1.25,
    }
    return params


def apply_params(plt, extra=None):
    p = get_params()
    if extra is not None:
        p.update(extra)
    plt.rcParams.update(p)


def adjust_ticks(ax, kx=0.25, ky=0.25, autoscale=True):
    def roundauto(x):
        dec = -int(np.floor(np.log10(x)))
        return np.round(x, dec)

    if kx is not None:
        ax.xaxis.set_major_locator(
            plticker.MultipleLocator(
                base=roundauto(np.ptp(ax.get_xticks()) * kx)))
    if ky is not None:
        ax.yaxis.set_major_locator(
            plticker.MultipleLocator(
                base=roundauto(np.ptp(ax.get_yticks()) * ky)))
    if autoscale:
        ax.autoscale()


def get_cycle(v, size=None):
    c = itertools.cycle(v)
    if size is None:
        return c
    if size == 0:
        return v
    return [next(c) for _ in range(size)]


def get_color_cycle(size=None):
    return get_cycle(colorscheme, size)


def get_marker_cycle(size=None):
    return get_cycle(('o', 's', 'v', '^', 'p', 'h', 'P'), size)


def get_line_cycle(size=None):
    return get_cycle(('-', '--', '-.', ':'), size)


def units(s):
    r"""
    Extracts space separated words and:
    * adds thin spaces if word is '/'
    * wraps by \mathrm{} if word consists of [a-zA-Z] and does not start with '\'
    Then removes spaces before '^'.

    Example:
    >>> units('Pa \cdot s / m ^3')
    '\mathrm{Pa} \cdot \mathrm{s} \,/\, \mathrm{m}^3'
    """
    words = s.split()
    words = [
        r"\,/\," if w == '/' else  #
        "\mathrm{" + w + "}" if not w[0] == '\\' and w.isalpha()  #
        else w for w in words
    ]
    res = ' '.join(words)
    res = res.replace(' ^', '^')
    return res


def bunits(s):
    return '$~[' + units(s) + ']$'


def get_parser(name, return_parser=False):
    parser = ArgumentParser()
    parser.add_argument('--conf',
                        default=name + ".conf.py",
                        help="path to configuration")
    parser.add_argument('--datadir',
                        default="",
                        help="path to data directory with data files"
                        ", defaults to directory containing configuration")
    parser.add_argument('--cachedir',
                        default="",
                        help="path to directory with plot data cache"
                        ", defaults to current directory")
    parser.add_argument('--update',
                        action="store_true",
                        help="process data and update cache")
    parser.add_argument('--outputdir',
                        default="",
                        help="path to output directory"
                        ", defaults to current directory")
    return parser


def parse_args(name=None, parser=None):
    if parser is None:
        parser = get_parser(name)
    args = parser.parse_args()
    if args.datadir == "":
        args.datadir = os.path.dirname(args.conf)
    return args


def execdict(path):
    conf = dict()
    with open(path) as f:
        exec(f.read(), conf)
    conf.pop('__builtins__')
    return Namespace(**conf)


def default_cache(args, name, arg0=False):
    return cache_to_file(os.path.join(args.cachedir, name + ".pickle"),
                         args.update, arg0)


def savefig(fig, path, **kwargs):
    print(path)
    fig.savefig(path, metadata={'CreationDate': None}, **kwargs)


def savelegend(fig, ax, path, detect_codes=False, **kwargs):
    """
    detect_codes: prepend labels with style codes detected from lines
    """
    figleg, axleg = plt.subplots()
    handles, labels = ax.get_legend_handles_labels()
    if detect_codes:
        labels = [
            code_to_str(line_to_code(h)) + ' ' + l
            for h, l in zip(handles, labels)
        ]
    legend = axleg.legend(handles,
                          labels,
                          loc='center',
                          frameon=False,
                          handlelength=2.68)
    axleg.axis('off')
    figleg.canvas.draw()
    bbox = legend.get_window_extent().transformed(
        fig.dpi_scale_trans.inverted())
    # FIXME workaround for many lines (>5), incorrect box size
    dy = bbox.y1 - bbox.y0
    bbox.y0 -= dy * 0.022
    bbox.y1 += dy * 0.022
    savefig(figleg, path, bbox_inches=bbox, **kwargs)

def get_styles():
    l = [''] + list(get_line_cycle(0))
    c = ['#000000'] + list(get_color_cycle(0))
    m = [''] + list(get_marker_cycle(0))
    return l, c, m


def get_style_codes():
    l, c, m = get_styles()
    return range(len(l)), range(len(c)), range(len(m))

def code_to_str(code):
    return ''.join([str(c) if c is not None else 'n' for c in code])

def code_to_style(code):
    """
    code: `str` or array-like
        Line code (line, color, marker).
        Examples: '123', [1, 2, 3]
    """
    if type(code) == str:
        code = list(map(int, code))
    l, c, m = get_styles()
    return l[code[0]], c[code[1]], m[code[2]]

def line_to_code(line):
    ll, cc, mm = get_styles()
    l = line.get_linestyle()
    c = line.get_color()
    m = line.get_marker()
    if l in [None, 'None']: l = ''
    if m in [None, 'None']: m = ''
    if c == 'k': c = '#000000'
    def index(tt, t):
        if t in tt:
            return tt.index(t)
        return None
    return index(ll, l), index(cc, c), index(mm, m)
