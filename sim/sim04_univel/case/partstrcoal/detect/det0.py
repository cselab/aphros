#!/usr/bin/env python

# Prints positions of the neck:
# zmin  , zmax

import numpy as np
import cv2
import os

import sys

av = sys.argv
f = av[1]  # input image

assert os.path.isfile(f)

im = cv2.imread(f)
gr = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)

u = gr.T.astype(float) / 255.
# image size
sx = u.shape[0]
sy = u.shape[1]
x = np.linspace(0, 1, sx)
y = np.linspace(0, 1, sy)

# threshold, inside the bubble if u < th
th = 0.5

v = u[:,-1]
ii = np.where(v < th)[0]

i0 = ii.min()
i1 = ii.max()

x0 = x[i0]
x1 = x[i1]

print("{:} {:}".format(x0, x1))
