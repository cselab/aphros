#!/usr/bin/env python

import argparse
import math
import os
import sys
import numpy as np

def NormalizeType(v):
    if isinstance(v, (float, np.float, np.float32, np.float64)):
        v = float(v)
    elif isinstance(v, (int, np.int, np.int32, np.int64)):
        v = int(v)
    elif isinstance(v, str):
        v = str(v)
    elif isinstance(v, list):
        v = list(map(float, v))
    else:
        assert False, "unknown type of value '{:}'".format(str(v))
    return v

class Var:
    def __init__(self):
        self.var = dict()
    def __getitem__(self, key):
        return self.var[key]
    def __setitem__(self, key, value):
        self.var[key] = value
    def Generate(self):
        t = []
        for k,v in self.var.items():
            v = NormalizeType(v)
            if isinstance(v, float):
                t.append("set double {:} {:}".format(k, v))
            elif isinstance(v, int):
                t.append("set int {:} {:}".format(k, v))
            elif isinstance(v, list):
                v = ' '.join(map(str, v))
                t.append("set vect {:} {:}".format(k, v))
            elif isinstance(v, str):
                if '\n' in v:
                    v = '"' + v + '"'
                t.append("set string {:} {:}".format(k, v))
            else:
                assert False, "unknown type of value '{:}'".format(str(v))
        return '\n'.join(t)

def VectToStr(v):
    return ' '.join(map(str, v))

class Geometry:
    def __init__(self):
        self.lines = []
    def __Append(self, line):
        self.lines.append(line)
    def __Prefix(self, kwargs):
        s = ''
        if kwargs.get("intersect", False):
            s += '&'
        if kwargs.get("invert", False):
            s += '-'
        return s
    def Box(self, center, halfsize, rotation_z=0, **kwargs):
        s = self.__Prefix(kwargs)
        s += "box {:} {:} {:}".format(
                VectToStr(center), VectToStr(halfsize), rotation_z)
        self.__Append(s)
    def Generate(self):
        return '\n'.join(self.lines)

class Bc:
    def __init__(self):
        self.lines = []
    def __Indent(self, text):
        lines = text.split('\n')
        return '\n'.join(["  " + line for line in lines])
    def __Append(self, line, geom):
        line = "{:} {{\n{:}\n}}".format(line, self.__Indent(geom.Generate()))
        self.lines.append(line)
    def Wall(self, geom, velocity):
        s = "wall {:}".format(VectToStr(velocity))
        self.__Append(s, geom)
    def SlipWall(self, geom):
        s = "slipwall"
        self.__Append(s, geom)
    def Generate(self):
        return '\n'.join(self.lines)


if __name__ == "__main__":
    var = Var()
    var["key0"] = 0
    var["key5"] = 0
    var["key5"] = 2
    var["key1"] = 0.
    var["key2"] = "str"
    var["key3"] = [0, 1, 2]
    print(var.Generate())
    print()

    geom = Geometry()
    geom.Box(center=[0,0,0], halfsize=[1,1,1], rotation_z=5)
    geom.Box(center=[0,2,2], halfsize=[1,1,1], rotation_z=5)
    print(geom.Generate())

    bc = Bc()
    bc.Wall(geom, velocity=[1,0,0])
    bc.SlipWall(geom)
    print(bc.Generate())
