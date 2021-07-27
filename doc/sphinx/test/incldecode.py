#!/usr/bin/env python3

import sys
import includecode

me = "test/includecode.py"
def usg():
    sys.stderr.write("%s [-s] name file\n" % me)
    sys.stderr.write("%s Simple /u/aphros/src/solver/simple.h\n" % me)
    sys.exit(1)

Get = includecode.GetFunc
if "-h" in sys.argv:
    usg()
if "-s" in sys.argv:
    Get = includecode.GetStruct
    sys.argv.remove("-s")
name = sys.argv[1]
filename = sys.argv[2]

with open(filename, 'r') as f:
    includecode.filename = filename
    f = f.read()
    print(Get(f, name))
