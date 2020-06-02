#!/usr/bin/env python3

# source: https://gist.github.com/pv/8036995

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.spatial import Voronoi, voronoi_plot_2d
import sys
import os

# Write uniform grid data
# u -- 2d or 3d array
# Format:
# <Nx> <Ny> <Nz>
# <u[0,0,0]> <u[1,0,0]> ...
def WritePlain(u, path):
    s = u.shape
    assert len(s) in [2, 3]
    if (len(s) == 2):
        u = u.reshape((s[0], s[1], 1))
    with open(path, 'w') as f:
        f.write("{:} {:} {:}\n".format(*u.shape))
        u = u.flatten()
        np.savetxt(f, u, newline='', fmt='%.16g ')

def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

# make up data points
'''
np.random.seed(1234)
nx = 10
ny = 10
gap = 0.1
xx = np.linspace(gap, 1 - gap, nx)
yy = np.linspace(gap, 1 - gap, ny)
points = np.meshgrid(xx, yy)
points = np.array(points)
points = points.reshape(2, nx * ny)
points = points.T
points += (np.random.rand(nx * ny, 2) - 0.5) * gap
'''

av = sys.argv

if len(av) < 2:
    sys.stderr.write('''./voronoi.py POINTS
Creates 2D Voronoi diagram from points.
POINTS: table with two columns (e.g. points0)
Output:
  %.dat: field in aphros plain data format (e.g. points0.dat)
  %.png: image of the diagram (e.g. points0.png)
''')
    exit(1)

points_path = av[1]
points = np.loadtxt(points_path)
print("Read {:} points from '{:}'".format(len(points), points_path))
vor = Voronoi(points)
basename = os.path.basename(points_path)

regions, vertices = voronoi_finite_polygons_2d(vor)

resx = 256
dpi = 100
fig = plt.figure(figsize=(resx/dpi,resx/dpi))
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_axis_off()
fig.add_axes(ax)

for i,region in enumerate(regions):
    polygon = vertices[region]
    c = "#{:02x}{:02x}{:02x}".format(i, i, i)
    ax.fill(*zip(*polygon), c=c, antialiased=False)

fig.canvas.draw()

data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
data = data[:,:,0]
data = np.flipud(data)
WritePlain(data, basename + ".dat")

voronoi_plot_2d(vor, ax, show_vertices=False, line_colors='#ffffff')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
fig.savefig(basename + ".png")
