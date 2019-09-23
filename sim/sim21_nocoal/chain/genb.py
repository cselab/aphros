#!/usr/bin/env python

import numpy as np

h = 0.075 * 2
r = 0.1
th = r / 3
dx = r + th


fn = "b.dat"

with open(fn, 'w') as f:
  z = dx
  while z < 1. - dx:
    y = dx
    while y < 1. - dx:
      x = dx
      while x < 1. - dx:
        f.write("ring {x} {y} {z} 1 0 0 {r} {th}\n".format(**locals()))
        x1 = x + dx
        y1 = y + dx
        z1 = z + dx
        f.write("ring {x} {y1} {z} 0 0 1 {r} {th}\n".format(**locals()))
        f.write("ring {x} {y} {z1} 0 1 0 {r} {th}\n".format(**locals()))

        f.write("ring {x1} {y1} {z} 0 1 0 {r} {th}\n".format(**locals()))
        f.write("ring {x1} {y} {z1} 0 0 1 {r} {th}\n".format(**locals()))

        f.write("ring {x1} {y1} {z1} 1 0 0 {r} {th}\n".format(**locals()))
        x += dx*2
      y += dx*2
    z += dx*2
