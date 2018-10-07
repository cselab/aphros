#!/usr/bin/env python

from glob import glob
import re
import ch.plot as chp
import numpy as np
import matplotlib as mp
import os

def ReadAll(f):
    return " ".join([l.strip() for l in open(f).readlines()])

# f: path to p*/trep_*
# Returns:
# n: number of cores
# t: execution time
def GetTrep(f):
    a = ReadAll(f)
    t = float(re.findall("all \[([0-9.]*)", a)[0])
    return t

def GetTotal(f):
    return int(re.findall("total = ([0-9]*)", a)[0])

# f: path to p*/out
# Returns:
# n: number of cores
# t: execution time
def GetTime0(f):
    fd = os.path.dirname(f)
    n = int(re.findall("p([0-9]*)", fd)[0])
    g = glob(os.path.join(fd, "trep_*"))
    tt = []
    for tr in g:
        t = GetTrep(tr)
        tt.append(GetTrep(tr))
    tt = np.array(tt)
    return n, np.median(tt)

def GetTime():
    ff = glob("p*/out")

    tt = dict()
    for f in ff:
        n,t = GetTime0(f)
        tt[n] = t

    return tt

# tt: dictionary tt {n:t}
# fo: output path
def WriteTime(tt, fo="time"):
    with open(fo, 'w') as o:
        o.write("n t\n")
        for n in sorted(tt):
            o.write("{:} {:}\n".format(n, tt[n]))

# Plot weak scaling
# tt: dictionary tt {n:t}
def PlotWeak(tt, po="weak.pdf"):
    # number of cores
    n = np.array(sorted(tt))
    # execution time
    t = np.array([tt[q] for q in n])
    # reference time
    ti = 0.*t + t[0]

    nd = n / 12

    fig, ax = chp.PlotInit()
    ax.plot(nd, ti / t, label='mfer', marker='.')
    ax.plot(nd, ti / ti, label='ideal')
    ax.set_xscale('log', basex=2)
    ax.set_xlabel("nodes")
    ax.set_ylabel("weak scaling efficiency")
    ax.minorticks_off()
    #ax.set_xticks([2 ** p for p in range(0, 20, 2)])
    #ax.set_xticklabels([2 ** p for p in range(0, 10, 2)])
    ax.legend()
    chp.PlotSave(fig, ax, po)

# Plot strong scaling
# tt: dictionary tt {n:t}
def PlotStrong(tt, po="strong.pdf", poe="strongeff.pdf"):
    # number of cores
    n = np.array(sorted(tt))
    # execution time
    t = np.array([tt[q] for q in n])
    # speedup
    s = t[0] / t * n[0]
    # ideal speedup
    si = n

    nd = n // 12

    fig, ax = chp.PlotInit()
    ax.plot(nd, s, label='mfer')
    ax.plot(nd, si, label='ideal')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("nodes")
    ax.set_ylabel("speedup")
    ax.minorticks_off()
    ax.set_xticks(n)
    ax.set_xticklabels(n)
    ax.set_yticks(n)
    ax.set_yticklabels(n)
    ax.legend()
    chp.PlotSave(fig, ax, po)

    fig, ax = chp.PlotInit()
    ax.plot(n, s / n, label='mfer')
    ax.plot(n, si / n, label='ideal')
    ax.set_xscale('log')
    ax.set_xlabel("cores")
    ax.set_ylabel("strong scaling efficiency")
    ax.minorticks_off()
    ax.set_xticks(n)
    ax.set_xticklabels(n)
    ax.legend()
    chp.PlotSave(fig, ax, poe)

tt = GetTime()
WriteTime(tt)
PlotWeak(tt)
