#!/usr/bin/env python

from ctypes import *

l = cdll.LoadLibrary("./libimp.so")


def callback(a, b, fcu):
    print(a, b, fcu)


c_callback = CFUNCTYPE(None, c_int, c_int, c_char_p)(callback)


l.CMain(c_callback)

