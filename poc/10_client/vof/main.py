#!/usr/bin/env python

from ctypes import *

l = cdll.LoadLibrary("./libimp.so")


def callback(a):
    print(a)


c_callback = CFUNCTYPE(None, c_int)(callback)


l.CMain(c_callback)

