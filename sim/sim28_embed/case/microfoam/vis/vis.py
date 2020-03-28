#!/usr/bin/env pvbatch


# state file generated using paraview version 5.8.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import paratools
import numpy as np

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


bcvtk = LegacyVTKReader(FileNames=paratools.FindFiles('bc.vtk'))
sm_0 = LegacyVTKReader(FileNames=paratools.FindFiles())
sources_ft, timearrays = paratools.ApplyForceTime([sm_0])
sm_0, = sources_ft

box = paratools.GetBox(bcvtk)
box = np.array(box)
boxc = (box[0] + box[1]) * 0.5
res = 1024
pixel = 1. / res
div = 16
boxsize = box[1] - box[0] + pixel * div
viewsize = [int(div * a * res + div - 1) // div for a in boxsize[:2]]

light1 = CreateLight()
light1.Position = [1.8562249839305878, 0.5329485405236483, 5.206477122922456]
light1.FocalPoint = [1.8562249839305878, 0.5329485405236483, 0.048663998022675514]
light1.Radius = 1.0

color_gray = [0.75] * 3

renderView1 = CreateView('RenderView')
renderView1.ViewSize = viewsize
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.CameraPosition = [boxc[0], boxc[1], 10]
renderView1.CameraFocalPoint = [boxc[0], boxc[1],  0]
renderView1.CameraParallelScale = 0.5 * max(viewsize) / (max(boxsize) * min(viewsize)) * 1.05
renderView1.CameraParallelProjection = 1
renderView1.Background = color_gray
renderView1.AdditionalLights = light1

SetActiveView(None)

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

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

slice1 = Slice(Input=bcvtk)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = boxc
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1.HyperTreeGridSlicer.Origin = [1.5499999523162842, 0.532812487334013, 0.04843749850988388]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from sm_0
sm_0Display = Show(sm_0, renderView1, 'GeometryRepresentation')
sm_0Display.Representation = 'Surface'
sm_0Display.AmbientColor = [0.0, 0.0, 0.0]
sm_0Display.ColorArrayName = ['POINTS', '']
sm_0Display.DiffuseColor = [0.0, 0.0, 0.0]
sm_0Display.Specular = 1.0
sm_0Display.SpecularColor = color_gray
sm_0Display.SpecularPower = 5.0
sm_0Display.Diffuse = 0.0
sm_0Display.Roughness = 0.2

slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
slice1Display.Representation = 'Wireframe'
slice1Display.AmbientColor = [0.0, 0.0, 0.0]
slice1Display.ColorArrayName = ['POINTS', '']
slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
slice1Display.LineWidth = 8.

paratools.SaveAnimation(renderView1, sources_ft, timearrays)
