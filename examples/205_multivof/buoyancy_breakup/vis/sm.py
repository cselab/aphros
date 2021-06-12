#!/usr/bin/env pvbatch

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

from glob import glob

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1080, 1080]
renderView1.InteractionMode = '2D'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 0.5, 0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.5, 0.5, 3.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.5
renderView1.CameraParallelProjection = 1
renderView1.Background = [1., 1., 1.]
renderView1.UseLight = 0

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

nx064 = LegacyVTKReader(FileNames=sorted(glob("nx064/sm_*.vtk")))
nx128 = LegacyVTKReader(FileNames=sorted(glob("nx128/sm_*.vtk")))
nx256 = LegacyVTKReader(FileNames=sorted(glob("nx256/sm_*.vtk")))

#d3 category10

vhex = [
    "1f77b4",
    "ff7f0e",
    "2ca02c",
    "d62728",
    "9467bd",
    "8c564b",
    "e377c2",
    "7f7f7f",
    "bcbd22",
    "17becf",
]


def rgb(h):
    return list(int(h[i:i + 2], 16) / 255. for i in (0, 2, 4))

lw = 8

nx064Display = Show(nx064, renderView1, 'GeometryRepresentation')
nx064Display.Representation = 'Wireframe'
nx064Display.AmbientColor = rgb(vhex[0])
nx064Display.ColorArrayName = ['POINTS', '']
nx064Display.DiffuseColor = rgb(vhex[0])
nx064Display.LineWidth = lw
nx064Display.Position = [0.0, 0.0, 0.0]

nx128Display = Show(nx128, renderView1, 'GeometryRepresentation')
nx128Display.Representation = 'Wireframe'
nx128Display.AmbientColor = rgb(vhex[1])
nx128Display.ColorArrayName = ['POINTS', '']
nx128Display.DiffuseColor = rgb(vhex[1])
nx128Display.LineWidth = lw
nx128Display.Position = [0.0, 0.0, 1.0]

nx256Display = Show(nx256, renderView1, 'GeometryRepresentation')
nx256Display.Representation = 'Wireframe'
nx256Display.AmbientColor = rgb(vhex[2])
nx256Display.ColorArrayName = ['POINTS', '']
nx256Display.DiffuseColor = rgb(vhex[2])
nx256Display.LineWidth = lw
nx256Display.Position = [0.0, 0.0, 2.0]

tk = GetTimeKeeper()


def save(t):
    print(t)
    tk.Time = t
    SaveScreenshot("a_{:04d}.png".format(t))


for t in range(0, 101, 10):
    save(t)
