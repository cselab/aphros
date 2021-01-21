#!/usr/bin/env python3

import trimesh

f = "a.stl"
m = trimesh.load_mesh(f, process=False)
msub = trimesh.remesh.subdivide_to_size(m.vertices, m.faces, max_edge=2)
msub = trimesh.Trimesh(*msub, process=False)
msub.export("a_remesh.off")
msub.export("a_remesh.stl")
