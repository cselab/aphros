#!/usr/bin/env python

# compute dk/dt : time derivative of kinetic energy

import numpy as np

# compute derivative, forward difference, extrapolate at last point
# u: function, 1d array
# t: intependent variable, 1d array
def GetDeriv(u, t):
  dt = np.roll(t, -1) - t
  dt[-1] = dt[-2]
  du = np.roll(u, -1) - u
  du[-1] = du[-2]
  return du / dt

# array from ascii column
# p: path
def ReadArray(p):
  u = np.loadtxt(p)
  return u.flatten()

# write array as ascii column
def WriteArray(u, p):
  np.savetxt(p, u)


d = "sc/"

k = ReadArray(d + "k")
t = ReadArray(d + "t")

vol = (np.pi * 2) ** 3

dk = GetDeriv(k, t)
dk = -dk / vol
WriteArray(dk, d + "dk")
