#!/usr/bin/env python3

import numpy as np
import aphros
import os


def V(x, y=None, z=None):
    if y is None:
        y = x
    if z is None:
        z = 0
    return np.array([x, y, z])


def Basis(angle_rad):
    '''
    Returns basis (eu,ev) produced by rotation of Cartasian basis by `angle_rad`
    '''
    eu = V(np.cos(angle_rad), np.sin(angle_rad))
    ev = V(-eu[1], eu[0])
    return eu, ev


def Dot(u, v):
    return (u * v).sum()

def Dist(u, v):
    return Dot(u - v, u - v) ** 0.5

def Cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def Intersection(pu, u, pv, v):
    '''
    Intersection point of two lines (pu,u) and (pv,v)
    specified by a point (pu) and direction (u)
    '''
    kv = Cross(u, pu - pv) / Cross(u, v)
    return pv + kv * v


def ClipBox(box_c, box_r, box_angle_deg, line_p, line_d):
    '''
    Clip box by line
    box_c, box_r: center and half-size
    line_p, line_d: line point and direction
    '''
    eu, ev = Basis(np.radians(box_angle_deg))
    inter = Intersection(box_c + ev * box_r[1], eu, line_p, line_d)
    corner = box_c - eu * box_r[0] + ev * box_r[1]
    dist = Dist(inter, corner)
    box_c = box_c + eu * dist * 0.5
    box_r = box_r + V(-dist * 0.5, 0)
    return box_c, box_r, box_angle_deg

def Valve(g, start, L, W, Ro, Ri, angle_rad, shift=V(0, 0)):
    '''
    Appends a Tesla valve to geometry.
    g: aphros.Geometry instance
    start: center of the left edge, (x,y)
    L: pipe length
    W: pipe width
    Ro: outer radius
    Ri: inner radius
    angle_rad: inclination angle (0 if horizontal) in radians
    shift: shift of circle position relative to position below the pipe end
    forward: pipe direction, True if forward
    Returns:
    end: center of the right edge, (x,y)
    '''

    eu, ev = Basis(angle_rad)

    Lh = L / 2
    Wh = W / 2
    angle_deg = np.degrees(angle_rad)


    end = start + eu * L

    inf = L * 10

    sgnx = 1 if eu[0] > 0 else -1
    sgny =  1 if eu[1] > 0 else -1

    shift = shift * V(sgnx, sgny)

    c = end + shift - V(Ro * sgnx, (Ro + Wh) * sgny)
    g.Cylinder(c, V(0, 0, 1), Ro, [-inf, inf])
    g.Cylinder(c, V(0, 0, 1), Ri, [-inf, inf], invert=True, intersect=True)

    eiu, eiv = Basis(-angle_rad)

    line = (start + ev * Wh, eu)

    b = (c - eiu * L + eiv * (Ro + Ri) / 2, V(L, (Ro - Ri) / 2), -angle_deg)
    b = ClipBox(*b, *line)
    g.Box(*b)

    b = (c - eiu * L - eiv * (Ro + Ri) / 2, V(L, (Ro - Ri) / 2), -angle_deg)
    b = ClipBox(*b, *line)
    g.Box(*b)

    b = (c - eiu * L, V(L, Ri), -angle_deg)
    b = ClipBox(*b, *line)
    g.Box(*b, invert=True, intersect=True)

    g.Box(start + eu * Lh, V(Lh, Wh), angle_deg)
    return end


class Default(aphros.Parameters):
    extent = 11.5
    angle_deg = 14
    W = 1.2
    Ro = 1.9
    Ri = Ro - W - 0.2
    L = 14
    cx = -1.0
    cy = 2.6
    shift = V(-2, -1.2)
    template = "tesla_template.conf"
    output = "tesla.conf"

path = "par.py"
if os.path.isfile(path):
    par = Default(path)
else:
    par = Default()

g = aphros.Geometry()

angle_rad = np.radians(par.angle_deg)

cx = par.cx
cy = par.cy
W = par.W
L = par.L
Ro = par.Ro
Ri = par.Ri
extent = par.extent

e = Valve(g, V(cx, cy), L, W=W, Ro=Ro, Ri=Ri, angle_rad=angle_rad, shift=par.shift)

e = Valve(g, V(extent-cx, extent - cy - 1.5), L, W=W, Ro=Ro, Ri=Ri, angle_rad=3.14-angle_rad, shift=par.shift)

body = g.Generate()

if par.template:
    with open(par.template) as f:
        t = f.read()
    t = t.replace("BODY", body)
    with open(par.output, 'w') as f:
        f.write(t)
else:
    g.GenerateFile(par.output)
