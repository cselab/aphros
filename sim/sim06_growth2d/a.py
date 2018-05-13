# state file generated using paraview version 5.4.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [963, 764]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 0.5, 0.06014088960364461]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.9407913785157281, -2.640372067346219, 1.0227406872687494]
renderView1.CameraFocalPoint = [0.5, 0.5000000000000002, 0.060140889603644374]
renderView1.CameraViewUp = [-0.028851832133593446, 0.28923928121341175, 0.9568219322244259]
renderView1.CameraParallelScale = 0.7088723439575985
renderView1.Background = [0.0, 0.0, 0.0]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
ppvd = PVDReader(FileName='/home/petr/s/mfer/ch/sim/sim06_growth2d/p.pvd')

# create a new 'Slice'
slice1 = Slice(Input=ppvd)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5, 0.5, 0.015625]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=slice1)

# create a new 'Warp By Scalar'
warpByScalar1 = WarpByScalar(Input=cellDatatoPointData1)
warpByScalar1.Scalars = ['POINTS', 'p']
warpByScalar1.ScaleFactor = 0.1

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'vf'
vfLUT = GetColorTransferFunction('vf')
vfLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'vf'
vfPWF = GetOpacityTransferFunction('vf')
vfPWF.ScalarRangeInitialized = 1

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.RGBPoints = [-1.6277899742126465, 0.231373, 0.298039, 0.752941, 3.147624969482422, 0.865003, 0.865003, 0.865003, 7.92303991317749, 0.705882, 0.0156863, 0.14902]
pLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')
pPWF.Points = [-1.6277899742126465, 0.0, 0.5, 0.0, 7.92303991317749, 1.0, 0.5, 0.0]
pPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'p']
slice1Display.LookupTable = pLUT
slice1Display.OSPRayScaleArray = 'p'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
slice1Display.GaussianRadius = 0.05
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# show data from warpByScalar1
warpByScalar1Display = Show(warpByScalar1, renderView1)
# trace defaults for the display properties.
warpByScalar1Display.Representation = 'Surface'
warpByScalar1Display.ColorArrayName = ['POINTS', 'vf']
warpByScalar1Display.LookupTable = vfLUT
warpByScalar1Display.OSPRayScaleArray = 'p'
warpByScalar1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByScalar1Display.SelectOrientationVectors = 'None'
warpByScalar1Display.ScaleFactor = 0.8382318824529649
warpByScalar1Display.SelectScaleArray = 'None'
warpByScalar1Display.GlyphType = 'Arrow'
warpByScalar1Display.GlyphTableIndexArray = 'None'
warpByScalar1Display.DataAxesGrid = 'GridAxesRepresentation'
warpByScalar1Display.PolarAxes = 'PolarAxesRepresentation'
warpByScalar1Display.GaussianRadius = 0.41911594122648244
warpByScalar1Display.SetScaleArray = ['POINTS', 'p']
warpByScalar1Display.ScaleTransferFunction = 'PiecewiseFunction'
warpByScalar1Display.OpacityArray = ['POINTS', 'p']
warpByScalar1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
warpByScalar1Display.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for pLUT in view renderView1
pLUTColorBar = GetScalarBar(pLUT, renderView1)
pLUTColorBar.Title = 'p'
pLUTColorBar.ComponentTitle = ''

# get color legend/bar for vfLUT in view renderView1
vfLUTColorBar = GetScalarBar(vfLUT, renderView1)
vfLUTColorBar.WindowLocation = 'UpperRightCorner'
vfLUTColorBar.Title = 'vf'
vfLUTColorBar.ComponentTitle = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(warpByScalar1)
# ----------------------------------------------------------------