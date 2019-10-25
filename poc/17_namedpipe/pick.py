#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
from numpy.random import rand
import threading
import sys


# Read uniform grid data
# lines
# Format:
# <nx> <ny> <nz>
# <u[z,y,x]> ...
# Return:
# array of shape (nx, ny, nz)
# None if file not found
def ReadPlain(ll):
    # shape x,y,z
    s = np.array(ll[0].split(), dtype=int)
    # shape z,y,x
    ss = tuple(reversed(s))
    # data flat
    u = np.array(ll[1].split(), dtype=float)
    # data z,y,x
    u = u.reshape(ss)
    return u


def Send(msg):
    try:
        with open('in', 'w') as f:
            f.write("{:}\n".format(msg))
    except IOError:
        print("can't open pipe in")

def Listen():
    while True:
        try:
            with open('out', 'r') as f:
                line = f.readline()
                if line.strip() == "field":
                    print(ReadPlain(f.readlines()))
                else:
                    print(f.readlines())
        except IOError:
            print("can't open pipe out")


def Print(msg, msg2=None):
    if msg2 is not None:
        msg = "{:}{:}".format(msg, msg2)
    print("send: " + msg)
    Send(msg)

def pick_simple():
    fig, ax = plt.subplots(1, 1)
    ax.text(0.0, 1.1, "exit", picker=True, bbox=dict(facecolor='red'), transform=ax.transAxes)
    ax.text(0.05, 1.1, "step", picker=True, bbox=dict(facecolor='green'), transform=ax.transAxes)
    ax.text(0.1, 1.1, "field", picker=True, bbox=dict(facecolor='yellow'), transform=ax.transAxes)
    line, = ax.plot(rand(100), 'o', picker=5)  # 5 points tolerance

    def onpick1(event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            Print("sum {:} {:}".format(xdata[ind][0], ydata[ind][0]))
        elif isinstance(event.artist, Rectangle):
            patch = event.artist
            Print('onpick1 patch:', patch.get_path())
        elif isinstance(event.artist, Text):
            text = event.artist
            Print(text.get_text())

    fig.canvas.mpl_connect('pick_event', onpick1)



if __name__ == '__main__':
    t1 = threading.Thread(target=Listen)
    t1.start()
    pick_simple()
    plt.show()
