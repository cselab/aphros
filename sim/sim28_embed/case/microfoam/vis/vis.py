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

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Light'
light1 = CreateLight()
light1.Position = [1.8562249839305878, 0.5329485405236483, 5.206477122922456]
light1.FocalPoint = [1.8562249839305878, 0.5329485405236483, 0.048663998022675514]
light1.Radius = 1.0

# get the material library
materialLibrary1 = GetMaterialLibrary()

# create light
# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2000, 700]
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.CameraPosition = [1.5524737638061836, 0.5325429962736342, 5.206477122922456]
renderView1.CameraFocalPoint = [1.5524737638061836, 0.5325429962736342, 0.048663998022675514]
renderView1.CameraParallelScale = 0.5509096758631732
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.7725490196078432, 0.7686274509803922, 0.7725490196078432]
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

# create a new 'Legacy VTK Reader'
bcvtk = LegacyVTKReader(FileNames=paratools.FindFiles('bc.vtk'))
sm_0 = LegacyVTKReader(FileNames=paratools.FindFiles())
sources_ft, timearrays = paratools.ApplyForceTime([sm_0])
sm_0, = sources_ft

box = paratools.GetBox(bcvtk)
boxc = [(box[0][i] + box[1][i]) * 0.5 for i in range(3)]
renderView1.CameraPosition[0] = boxc[0]
renderView1.CameraPosition[1] = boxc[1]
renderView1.CameraFocalPoint[0] = boxc[0]
renderView1.CameraFocalPoint[1] = boxc[1]

# create a new 'Plane'
plane1 = Plane()
plane1.Origin = [0.0, -1.0, 0.0]
plane1.Point1 = [4.0, -1.0, 0.0]
plane1.Point2 = [0.0, 2.0, 0.0]

# create a new 'Slice'
slice1 = Slice(Input=bcvtk)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = boxc
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [1.5499999523162842, 0.532812487334013, 0.04843749850988388]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from sm_0
sm_0Display = Show(sm_0, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
sm_0Display.Representation = 'Surface'
sm_0Display.AmbientColor = [0.0, 0.0, 0.0]
sm_0Display.ColorArrayName = ['POINTS', '']
sm_0Display.DiffuseColor = [0.0, 0.0, 0.0]
sm_0Display.PointSize = 30.0
sm_0Display.LineWidth = 3.0
sm_0Display.RenderPointsAsSpheres = 1
sm_0Display.Specular = 1.0
sm_0Display.SpecularColor = [0.7607843137254902, 0.7568627450980392, 0.7607843137254902]
sm_0Display.SpecularPower = 5.0
sm_0Display.Diffuse = 0.0
sm_0Display.Roughness = 0.2
sm_0Display.OSPRayScaleArray = 'nn'
sm_0Display.OSPRayScaleFunction = 'PiecewiseFunction'
sm_0Display.OSPRayMaterial = 'thin glass'
sm_0Display.SelectOrientationVectors = 'nn'
sm_0Display.ScaleFactor = 0.2475440561771393
sm_0Display.SelectScaleArray = 'c'
sm_0Display.GlyphType = 'Arrow'
sm_0Display.GlyphTableIndexArray = 'c'
sm_0Display.GaussianRadius = 0.012377202808856964
sm_0Display.SetScaleArray = ['POINTS', 'nn']
sm_0Display.ScaleTransferFunction = 'PiecewiseFunction'
sm_0Display.OpacityArray = ['POINTS', 'nn']
sm_0Display.OpacityTransferFunction = 'PiecewiseFunction'
sm_0Display.DataAxesGrid = 'GridAxesRepresentation'
sm_0Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sm_0Display.ScaleTransferFunction.Points = [-0.9931906461715698, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sm_0Display.OpacityTransferFunction.Points = [-0.9931906461715698, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from plane1
plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plane1Display.Representation = 'Surface'
plane1Display.AmbientColor = [0.7725490196078432, 0.7686274509803922, 0.7725490196078432]
plane1Display.ColorArrayName = [None, '']
plane1Display.DiffuseColor = [0.7725490196078432, 0.7686274509803922, 0.7725490196078432]
plane1Display.PointSize = 30.0
plane1Display.LineWidth = 3.0
plane1Display.RenderPointsAsSpheres = 1
plane1Display.OSPRayScaleArray = 'Normals'
plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane1Display.SelectOrientationVectors = 'None'
plane1Display.ScaleFactor = 0.30000000000000004
plane1Display.SelectScaleArray = 'None'
plane1Display.GlyphType = 'Arrow'
plane1Display.GlyphTableIndexArray = 'None'
plane1Display.GaussianRadius = 0.015
plane1Display.SetScaleArray = ['POINTS', 'Normals']
plane1Display.ScaleTransferFunction = 'PiecewiseFunction'
plane1Display.OpacityArray = ['POINTS', 'Normals']
plane1Display.OpacityTransferFunction = 'PiecewiseFunction'
plane1Display.DataAxesGrid = 'GridAxesRepresentation'
plane1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plane1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plane1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Wireframe'
slice1Display.AmbientColor = [0.0, 0.0, 0.0]
slice1Display.ColorArrayName = ['POINTS', '']
slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
slice1Display.PointSize = 30.0
slice1Display.LineWidth = 3.0
slice1Display.RenderPointsAsSpheres = 1
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.3099999904632569
slice1Display.SelectScaleArray = 'group'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'group'
slice1Display.GaussianRadius = 0.015499999523162842
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

paratools.SaveAnimation(renderView1, sources_ft, timearrays)
