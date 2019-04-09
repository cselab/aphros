#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from skimage import measure
from skimage.draw import ellipsoid

import h5py

fn="vf_0040.h5"
fn="vf_0000.h5"

f = h5py.File(fn)
u = np.array(f["data"])
f.close()

# Generate a level set about zero of two identical ellipsoids in 3D
#ellip_base = ellipsoid(6, 10, 16, levelset=True)
#ellip_double = np.concatenate((ellip_base[:-1, ...],
#                               ellip_base[2:, ...]), axis=0)

u = u[:,:,:,0]
s = 1
u = u[::s,::s,::s]
u = np.concatenate((u[:,:,:], u[:,::-1,:]), axis=1)
u = np.concatenate((u[:,:,:], u[:,:,::-1]), axis=2)

nz,ny,nx = u.shape
u = u[:nz//3, ny//4:ny//4*3, nx//5:nx*4//5]
hx = 1. / nz

u = np.transpose(u)

nz,ny,nx = u.shape

print(u.shape)
#exit()

# Use marching cubes to obtain the surface mesh of these ellipsoids
verts, faces, normals, values = \
    measure.marching_cubes_lewiner(u, 0.5)

area = measure.mesh_surface_area(verts, faces)
print(area * hx ** 2)

# Display resulting triangular mesh using Matplotlib. This can also be done
# with mayavi (see skimage.measure.marching_cubes_lewiner docstring).
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Fancy indexing: `verts[faces]` to generate a collection of triangles
mesh = Poly3DCollection(verts[faces])
mesh.set_edgecolor('k')
ax.add_collection3d(mesh)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax.azim = 80
ax.elev = 10

ax.set_aspect("equal")

nn = max(nx, ny, nz)
ax.set_xlim(0, nn)  # a = 6 (times two for 2nd ellipsoid)
ax.set_ylim(0, nn)  # b = 10
ax.set_zlim(0, nn)  # c = 16

plt.tight_layout()
plt.savefig("a.pdf")
