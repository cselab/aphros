#!/usr/bin/env python

import h5py
import numpy as np
#import matplotlib.pyplot as plt
import sys
import os

av = sys.argv

if len(av) not in [2] or av[1] == '-h':
  print("""usage: ./{:} HDF
""".format(os.path.basename(av[0])))
  exit(1)

# filename
fn = av[1]

f = h5py.File(fn, 'r')

u = f['data']

u = np.array(u)

f.close()

nx,ny,nz = u.shape[0:3]

u = u[:,:,:,0]

print(u[:,0,:].mean() - u[:,-1,:].mean())

#plt.imshow(u[:,0,:])
#plt.show()
