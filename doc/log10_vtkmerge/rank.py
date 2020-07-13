#!/usr/bin/env python

import numpy as np
import aphros
from collections import defaultdict

path = "sm_0000.vtk"
points,poly,fields = aphros.ReadVtkPoly(path)

ranks = defaultdict(int)

for p in poly:
    for i in p:
        ranks[i] += 1

counts = defaultdict(int)

for i,r in ranks.items():
    counts[r] += 1

for r in sorted(counts):
    print(r, counts[r])

with open("a.csv", 'w') as f:
    def W(l):
        f.write(str(l) + '\n')
    W('x,y,z,rank')
    for i in ranks:
        W('{:},{:},{:},{:}'.format(*points[i], ranks[i]))
