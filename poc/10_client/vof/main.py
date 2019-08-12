#!/usr/bin/env python

from ctypes import *
import numpy as np

l = cdll.LoadLibrary("./libimp.so")


def callback(nx, ny, nz, fcu):
    print(nx, ny, nz, fcu)
    s = (nz, ny, nx)
    u = np.ctypeslib.as_array(fcu,shape=s)
    u[:,3,3] += 1
    print(u)


c_callback = CFUNCTYPE(
        None, c_int, c_int, c_int, POINTER(c_double))(callback)


l.CMain(c_callback)

