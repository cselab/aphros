#!/usr/bin/env python

from ctypes import *
import numpy as np
import matplotlib.pyplot as plt
import sys

l = cdll.LoadLibrary("./libimp.so")

def callback(nx, ny, nz, fcu):
    print(nx, ny, nz, fcu)
    s = (nz, ny, nx)
    u = np.ctypeslib.as_array(fcu,shape=s)
    u = u.copy()
    print(u.mean())
    print(np.median(u))
    plt.imshow(u[3,:,:])
    plt.show()

TF = CFUNCTYPE(None, c_int, c_int, c_int, POINTER(c_double))
c_callback = TF(callback)


LP_c_char = POINTER(c_char)
LP_LP_c_char = POINTER(LP_c_char)

argc = len(sys.argv)
argv = (LP_c_char * (argc + 1))()
for i, arg in enumerate(sys.argv):
    e = arg.encode('utf-8')
    argv[i] = create_string_buffer(e)

l.CMain.argtypes = (c_int, LP_LP_c_char,TF)
l.CMain(argc, argv, c_callback)
