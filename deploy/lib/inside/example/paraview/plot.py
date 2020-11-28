#!/usr/bin/env pvbatch

from paraview.simple import *

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1080, 1080]
renderView1.OrientationAxesVisibility = 1
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]

renderView1.CameraPosition = \
[3.009763730471278, 1.4254855945586373, 1.5307173395032143]
renderView1.CameraFocalPoint =\
[0.12253547594266075, 0.33964250657439565, 0.23417955862837347]
renderView1.CameraViewUp = \
[-0.35273228579947324, -0.16242352949554179, 0.9215196859649226]
renderView1.Background = [1.0, 1.0, 1.0]

surf = LegacyVTKReader(FileNames=['bc.vtk'])

surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = "Surface With Edges"
surfDisplay.ColorArrayName = [None, '']
surfDisplay.DiffuseColor = [1.0, 1.0, 1.0]

fn = "a.png"
print(fn)
SaveScreenshot(fn)
