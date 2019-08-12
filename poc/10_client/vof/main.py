#!/usr/bin/env python

from ctypes import *
import numpy as np

l = cdll.LoadLibrary("./libimp.so")


def callback(a, b, fcu):
    print(a, b, fcu)
    u = np.ctypeslib.as_array(fcu,shape=(3,))
    print(u)


c_callback = CFUNCTYPE(
        None, c_int, c_int, POINTER(c_double))(callback)


l.CMain(c_callback)

