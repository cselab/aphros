#!/usr/bin/env python

import h5py
import numpy as np
import sys
import os

av = sys.argv

if len(av) not in [2] or av[1] == '-h':
  print("""usage: ./{:} HDF
Computes mean values over slices z=0 and z=nz-1.
HDF: data of shape (nz,ny,nz,1)
""".format(os.path.basename(av[0])))
  exit(1)

# filename
fn = av[1]

f = h5py.File(fn, 'r')

u = f['data']

u = np.array(u)

f.close()

u = u[:,:,:,0]

print("{:} {:}".format(u[0,:,:].mean(), u[-1,:,:].mean()))
