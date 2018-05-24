# state file generated using paraview version 5.5.0-RC3

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.0-RC3

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1573, 901]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 0.5, 0.0078125]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.467998152246415, 0.511568248550334, 2.73986330756888]
renderView1.CameraFocalPoint = [0.467998152246415, 0.511568248550334, 0.0078125]
renderView1.CameraParallelScale = 0.170440119823333
renderView1.Background = [0.32, 0.34, 0.43]

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Structured Grid Reader'
p48vts = XMLStructuredGridReader(FileName=['p.48.vts'])
p48vts.CellArrayStatus = ['vx', 'vy', 'vz', 'p', 'vf']

# create a new 'CSV Reader'
partitcsv = CSVReader(FileName=['partit.0.csv', 'partit.1.csv', 'partit.2.csv', 'partit.3.csv', 'partit.4.csv', 'partit.5.csv', 'partit.6.csv', 'partit.7.csv', 'partit.8.csv', 'partit.9.csv', 'partit.10.csv', 'partit.11.csv', 'partit.12.csv', 'partit.13.csv', 'partit.14.csv', 'partit.15.csv', 'partit.16.csv', 'partit.17.csv', 'partit.18.csv', 'partit.19.csv', 'partit.20.csv', 'partit.21.csv', 'partit.22.csv', 'partit.23.csv', 'partit.24.csv', 'partit.25.csv', 'partit.26.csv', 'partit.27.csv', 'partit.28.csv', 'partit.29.csv', 'partit.30.csv', 'partit.31.csv', 'partit.32.csv', 'partit.33.csv', 'partit.34.csv', 'partit.35.csv', 'partit.36.csv', 'partit.37.csv', 'partit.38.csv', 'partit.39.csv', 'partit.40.csv', 'partit.41.csv', 'partit.42.csv', 'partit.43.csv', 'partit.44.csv', 'partit.45.csv', 'partit.46.csv', 'partit.47.csv', 'partit.48.csv', 'partit.49.csv', 'partit.50.csv', 'partit.51.csv', 'partit.52.csv', 'partit.53.csv', 'partit.54.csv', 'partit.55.csv', 'partit.56.csv', 'partit.57.csv', 'partit.58.csv', 'partit.59.csv', 'partit.60.csv', 'partit.61.csv', 'partit.62.csv', 'partit.63.csv', 'partit.64.csv', 'partit.65.csv', 'partit.66.csv', 'partit.67.csv', 'partit.68.csv', 'partit.69.csv', 'partit.70.csv', 'partit.71.csv', 'partit.72.csv', 'partit.73.csv', 'partit.74.csv', 'partit.75.csv', 'partit.76.csv', 'partit.77.csv', 'partit.78.csv', 'partit.79.csv', 'partit.80.csv', 'partit.81.csv', 'partit.82.csv', 'partit.83.csv', 'partit.84.csv', 'partit.85.csv', 'partit.86.csv', 'partit.87.csv', 'partit.88.csv', 'partit.89.csv', 'partit.90.csv', 'partit.91.csv', 'partit.92.csv', 'partit.93.csv', 'partit.94.csv', 'partit.95.csv', 'partit.96.csv', 'partit.97.csv', 'partit.98.csv', 'partit.99.csv'])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=partitcsv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'
tableToPoints1.KeepAllDataArrays = 1

# create a new 'Slice'
slice1 = Slice(Input=p48vts)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5, 0.5, 0.0078125]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from tableToPoints1
tableToPoints1Display = Show(tableToPoints1, renderView1)

# get color transfer function/color map for 'c'
cLUT = GetColorTransferFunction('c')
cLUT.EnableOpacityMapping = 1
cLUT.RGBPoints = [0.0, 1.0, 0.0, 0.0, 2.500005, 1.0, 0.0, 1.0, 5.0000025, 0.0, 0.0, 1.0, 7.5, 0.0, 1.0, 1.0, 9.999975, 0.0, 1.0, 0.0, 12.49995, 1.0, 1.0, 0.0, 15.0, 1.0, 0.0, 0.0]
cLUT.ColorSpace = 'RGB'
cLUT.NanColor = [1.0, 0.0, 0.0]
cLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
tableToPoints1Display.Representation = 'Points'
tableToPoints1Display.ColorArrayName = ['POINTS', 'c']
tableToPoints1Display.LookupTable = cLUT
tableToPoints1Display.PointSize = 3.0
tableToPoints1Display.OSPRayScaleArray = 'c'
tableToPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display.SelectOrientationVectors = 'None'
tableToPoints1Display.ScaleFactor = 0.059375
tableToPoints1Display.SelectScaleArray = 'None'
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.GlyphTableIndexArray = 'None'
tableToPoints1Display.GaussianRadius = 0.0296875
tableToPoints1Display.CustomShader = ''
tableToPoints1Display.SetScaleArray = ['POINTS', 'c']
tableToPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.OpacityArray = ['POINTS', 'c']
tableToPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display.SelectionCellLabelFontFile = ''
tableToPoints1Display.SelectionPointLabelFontFile = ''
tableToPoints1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tableToPoints1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 15.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tableToPoints1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 15.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
tableToPoints1Display.DataAxesGrid.XTitleFontFile = ''
tableToPoints1Display.DataAxesGrid.YTitleFontFile = ''
tableToPoints1Display.DataAxesGrid.ZTitleFontFile = ''
tableToPoints1Display.DataAxesGrid.XLabelFontFile = ''
tableToPoints1Display.DataAxesGrid.YLabelFontFile = ''
tableToPoints1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
tableToPoints1Display.PolarAxes.PolarAxisTitleFontFile = ''
tableToPoints1Display.PolarAxes.PolarAxisLabelFontFile = ''
tableToPoints1Display.PolarAxes.LastRadialAxisTextFontFile = ''
tableToPoints1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from slice1
slice1Display = Show(slice1, renderView1)

# get color transfer function/color map for 'vf'
vfLUT = GetColorTransferFunction('vf')
vfLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.6, 0.6, 0.6]
vfLUT.ColorSpace = 'RGB'
vfLUT.AboveRangeColor = [1.0, 1.0, 1.0]
vfLUT.NanColor = [1.0, 0.0, 0.0]
vfLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'vf']
slice1Display.LookupTable = vfLUT
slice1Display.OSPRayScaleArray = 'p'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.05
slice1Display.CustomShader = ''
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.SelectionCellLabelFontFile = ''
slice1Display.SelectionPointLabelFontFile = ''
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
slice1Display.DataAxesGrid.XTitleFontFile = ''
slice1Display.DataAxesGrid.YTitleFontFile = ''
slice1Display.DataAxesGrid.ZTitleFontFile = ''
slice1Display.DataAxesGrid.XLabelFontFile = ''
slice1Display.DataAxesGrid.YLabelFontFile = ''
slice1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
slice1Display.PolarAxes.PolarAxisTitleFontFile = ''
slice1Display.PolarAxes.PolarAxisLabelFontFile = ''
slice1Display.PolarAxes.LastRadialAxisTextFontFile = ''
slice1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from p48vts
p48vtsDisplay = Show(p48vts, renderView1)

# trace defaults for the display properties.
p48vtsDisplay.Representation = 'Outline'
p48vtsDisplay.ColorArrayName = [None, '']
p48vtsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
p48vtsDisplay.SelectOrientationVectors = 'None'
p48vtsDisplay.ScaleFactor = 0.1
p48vtsDisplay.SelectScaleArray = 'None'
p48vtsDisplay.GlyphType = 'Arrow'
p48vtsDisplay.GlyphTableIndexArray = 'None'
p48vtsDisplay.GaussianRadius = 0.005
p48vtsDisplay.SetScaleArray = [None, '']
p48vtsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
p48vtsDisplay.OpacityArray = [None, '']
p48vtsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
p48vtsDisplay.DataAxesGrid = 'GridAxesRepresentation'
p48vtsDisplay.SelectionCellLabelFontFile = ''
p48vtsDisplay.SelectionPointLabelFontFile = ''
p48vtsDisplay.PolarAxes = 'PolarAxesRepresentation'
p48vtsDisplay.ScalarOpacityUnitDistance = 0.0883937422803018

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
p48vtsDisplay.DataAxesGrid.XTitleFontFile = ''
p48vtsDisplay.DataAxesGrid.YTitleFontFile = ''
p48vtsDisplay.DataAxesGrid.ZTitleFontFile = ''
p48vtsDisplay.DataAxesGrid.XLabelFontFile = ''
p48vtsDisplay.DataAxesGrid.YLabelFontFile = ''
p48vtsDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
p48vtsDisplay.PolarAxes.PolarAxisTitleFontFile = ''
p48vtsDisplay.PolarAxes.PolarAxisLabelFontFile = ''
p48vtsDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
p48vtsDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for vfLUT in view renderView1
vfLUTColorBar = GetScalarBar(vfLUT, renderView1)
vfLUTColorBar.Title = 'vf'
vfLUTColorBar.ComponentTitle = ''
vfLUTColorBar.TitleFontFile = ''
vfLUTColorBar.LabelFontFile = ''

# set color bar visibility
vfLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'c'
cPWF = GetOpacityTransferFunction('c')
cPWF.Points = [0.0, 0.0, 0.5, 0.0, 8.5248441696167, 0.0, 0.5, 0.0, 8.5248441696167, 1.0, 0.5, 0.0, 10.8074531555176, 1.0, 0.5, 0.0, 10.8074531555176, 0.0, 0.5, 0.0, 15.0, 0.0, 0.5, 0.0]
cPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'vf'
vfPWF = GetOpacityTransferFunction('vf')
vfPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(slice1)
# ----------------------------------------------------------------
