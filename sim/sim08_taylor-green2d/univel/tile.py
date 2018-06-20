#!/usr/bin/env python

import sys
import subprocess as sp

usg = '''tile list of movies
./tile.py nx ny [nx*ny paths, row major] out
'''

av = sys.argv
m = 1

if len(av) <= 1:
    print(usg)
    exit(1)

nx = int(av[m]); m += 1
ny = int(av[m]); m += 1

n = nx * ny

print(n)
vs = []
for k in range(n):
    vs.append(av[m]); m += 1

out = av[m]; m += 1

rx = 1920
ry = 1080

sx = rx // nx
sy = ry // ny

sx = min(sx, sy)
sy = sx

rx = sx * nx
ry = sy * ny

o = "ffmpeg"
for k in range(n):
    o += " -i {:}".format(vs[k])

o += " -filter_complex 'nullsrc=size={:}x{:} [base]; ".format(rx, ry)
for k in range(n):
    o += "[{:}:v] setpts=PTS-STARTPTS, scale={:}x{:} [b{:}]; ".format(k, sx, sy, k)

pr = "base"
for j in range(ny):
    for i in range(nx):
        k = j * nx + i
        b = "b{:}".format(k)
        t = "t{:}".format(k)
        o += "[{:}][{:}] overlay=shortest=1:x={:}:y={:}".format(
                pr, b, sx * i, sy * j)
        if k + 1 < n:
            o += " [{:}]; ".format(t)
        pr = t
o += "' -c:v libx264 {:}".format(out)

print(o)

sp.check_call(o, shell=True)




