#!/usr/bin/env python

from ctypes import *

l = cdll.LoadLibrary("./libimp.so")


def callback():
    print(1)


c_callback = CFUNCTYPE(None)(callback)


l.CMain(c_callback)

