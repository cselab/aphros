#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
from numpy.random import rand
import threading
import sys


field = None
fig = None
ax = None

# Read uniform grid data
# nx ny nz
# u[z,y,x] ...
def ReadPlain(ll):
    shape = np.array(ll[0].split(), dtype=int)
    shape = tuple(reversed(shape))
    u = np.array(ll[1].split(), dtype=float)
    u = u.reshape(shape)
    return u


def Send(msg):
    try:
        with open('in', 'w') as f:
            f.write("{:}\n".format(msg))
    except IOError:
        print("can't open pipe in")


def Listen():
    global field
    while True:
        try:
            with open('out', 'r') as f:
                line = f.readline()
                if line.strip() == "field":
                    field = ReadPlain(f.readlines())
                    ax.imshow(field[0,:,:])
                    fig.canvas.draw()
                else:
                    sys.stdout.write(f.readline())
        except IOError:
            print("can't open pipe out")


def Print(msg, msg2=None):
    if msg2 is not None:
        msg = "{:}{:}".format(msg, msg2)
    print("send: " + msg)
    Send(msg)

def pick_simple():
    global field
    global fig, ax
    fig, ax = plt.subplots(1, 1)
    ax.text(0.1, 0.9, "exit", picker=True, bbox=dict(facecolor='red'), transform=fig.transFigure)
    ax.text(0.15, 0.9, "step", picker=True, bbox=dict(facecolor='green'), transform=fig.transFigure)
    ax.text(0.2, 0.9, "field", picker=True, bbox=dict(facecolor='yellow'), transform=fig.transFigure)
    ax.text(0.25, 0.9, "step+field", picker=True, bbox=dict(facecolor='yellow'), transform=fig.transFigure)
    #line, = ax.plot(rand(100), 'o', picker=5)  # 5 points tolerance

    def onpick(event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            Print("sum {:} {:}".format(xdata[ind][0], ydata[ind][0]))
        elif isinstance(event.artist, Rectangle):
            patch = event.artist
            Print('onpick patch:', patch.get_path())
        elif isinstance(event.artist, Text):
            text = event.artist.get_text()
            if text == "step+field":
                Print("step")
                Print("field")
            else:
                Print(text)

    fig.canvas.mpl_connect('pick_event', onpick)



if __name__ == '__main__':
    t1 = threading.Thread(target=Listen)
    t1.start()
    pick_simple()
    plt.show()
