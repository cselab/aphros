# state file generated using paraview version 5.5.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Bar Chart View'
barChartView1 = CreateView('XYBarChartView')
barChartView1.ViewSize = [989, 1205]
barChartView1.ChartTitleFontFile = ''
barChartView1.LeftAxisTitleFontFile = ''
barChartView1.LeftAxisRangeMaximum = 32000.0
barChartView1.LeftAxisLabelFontFile = ''
barChartView1.BottomAxisTitleFontFile = ''
barChartView1.BottomAxisRangeMinimum = -15.0
barChartView1.BottomAxisRangeMaximum = 30.0
barChartView1.BottomAxisLabelFontFile = ''
barChartView1.RightAxisRangeMaximum = 6.66
barChartView1.RightAxisLabelFontFile = ''
barChartView1.TopAxisTitleFontFile = ''
barChartView1.TopAxisRangeMaximum = 6.66
barChartView1.TopAxisLabelFontFile = ''

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [990, 1205]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.9835446619468228, 0.8477150097571433, 0.3384675507214473]
renderView1.StereoType = 0
renderView1.CameraPosition = [3.071169240119643, 1.9297821164605595, 0.04016103178915791]
renderView1.CameraFocalPoint = [-0.5038006819805507, 0.10206901912342722, 0.5610656794830171]
renderView1.CameraViewUp = [0.29781795128211047, -0.7506224930483468, -0.5898053414678532]
renderView1.CameraParallelScale = 0.8660254037844386
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.Shadows = 1
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'CSV Reader'
partit0csv = CSVReader(FileName=['/home/kpetr/s/ch_partstr/curvbub/partit.0.csv'])

# create a new 'Legacy VTK Reader'
s0vtk = LegacyVTKReader(FileNames=['/home/kpetr/s/ch_partstr/curvbub/s.0.vtk'])

# create a new 'CSV Reader'
acsv = CSVReader(FileName=['/home/kpetr/s/ch_partstr/curvbub/a.csv'])

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=acsv)
tableToPoints2.XColumn = ''
tableToPoints2.YColumn = ''
tableToPoints2.ZColumn = ''

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=partit0csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Histogram'
histogram1 = Histogram(Input=tableToPoints1)
histogram1.SelectInputArray = ['POINTS', 'k']
histogram1.BinCount = 50
histogram1.CustomBinRanges = [11.676985368207644, 17.96547120173813]

# ----------------------------------------------------------------
# setup the visualization in view 'barChartView1'
# ----------------------------------------------------------------

# show data from histogram1
histogram1Display = Show(histogram1, barChartView1)

# trace defaults for the display properties.
histogram1Display.CompositeDataSetIndex = [0]
histogram1Display.AttributeType = 'Row Data'
histogram1Display.UseIndexForXAxis = 0
histogram1Display.XArrayName = 'bin_extents'
histogram1Display.SeriesVisibility = ['bin_values']
histogram1Display.SeriesLabel = ['bin_extents', 'bin_extents', 'bin_values', 'bin_values']
histogram1Display.SeriesColor = ['bin_extents', '0', '0', '0', 'bin_values', '0.89', '0.1', '0.11']
histogram1Display.SeriesPlotCorner = ['bin_extents', '0', 'bin_values', '0']
histogram1Display.SeriesLabelPrefix = ''

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from s0vtk
s0vtkDisplay = Show(s0vtk, renderView1)

# trace defaults for the display properties.
s0vtkDisplay.Representation = 'Surface'
s0vtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
s0vtkDisplay.ColorArrayName = [None, '']
s0vtkDisplay.Opacity = 0.55
s0vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s0vtkDisplay.SelectOrientationVectors = 'None'
s0vtkDisplay.ScaleFactor = 0.1
s0vtkDisplay.SelectScaleArray = 'None'
s0vtkDisplay.GlyphType = 'Arrow'
s0vtkDisplay.GlyphTableIndexArray = 'None'
s0vtkDisplay.GaussianRadius = 0.005
s0vtkDisplay.SetScaleArray = [None, '']
s0vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
s0vtkDisplay.OpacityArray = [None, '']
s0vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
s0vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
s0vtkDisplay.SelectionCellLabelFontFile = ''
s0vtkDisplay.SelectionPointLabelFontFile = ''
s0vtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s0vtkDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
s0vtkDisplay.DataAxesGrid.XTitleFontFile = ''
s0vtkDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
s0vtkDisplay.DataAxesGrid.YTitleFontFile = ''
s0vtkDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
s0vtkDisplay.DataAxesGrid.ZTitleFontFile = ''
s0vtkDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
s0vtkDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
s0vtkDisplay.DataAxesGrid.XLabelFontFile = ''
s0vtkDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
s0vtkDisplay.DataAxesGrid.YLabelFontFile = ''
s0vtkDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
s0vtkDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s0vtkDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
s0vtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
s0vtkDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
s0vtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
s0vtkDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
s0vtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
s0vtkDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
s0vtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from tableToPoints1
tableToPoints1Display = Show(tableToPoints1, renderView1)

# get color transfer function/color map for 'k'
kLUT = GetColorTransferFunction('k')
kLUT.AutomaticRescaleRangeMode = 'Never'
kLUT.EnableOpacityMapping = 1
kLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 14.0, 0.865003, 0.865003, 0.865003, 28.0, 0.705882, 0.0156863, 0.14902]
kLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
tableToPoints1Display.Representation = 'Points'
tableToPoints1Display.AmbientColor = [0.0, 0.0, 0.0]
tableToPoints1Display.ColorArrayName = ['POINTS', 'k']
tableToPoints1Display.LookupTable = kLUT
tableToPoints1Display.PointSize = 5.0
tableToPoints1Display.RenderPointsAsSpheres = 1
tableToPoints1Display.OSPRayScaleArray = 'c'
tableToPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display.SelectOrientationVectors = 'None'
tableToPoints1Display.ScaleFactor = 0.1073496940328927
tableToPoints1Display.SelectScaleArray = 'None'
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.GlyphTableIndexArray = 'None'
tableToPoints1Display.GaussianRadius = 0.005367484701644635
tableToPoints1Display.SetScaleArray = ['POINTS', 'c']
tableToPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.OpacityArray = ['POINTS', 'c']
tableToPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display.SelectionCellLabelFontFile = ''
tableToPoints1Display.SelectionPointLabelFontFile = ''
tableToPoints1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tableToPoints1Display.ScaleTransferFunction.Points = [14015.0, 0.0, 0.5, 0.0, 31027025.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tableToPoints1Display.OpacityTransferFunction.Points = [14015.0, 0.0, 0.5, 0.0, 31027025.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
tableToPoints1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
tableToPoints1Display.DataAxesGrid.XTitleFontFile = ''
tableToPoints1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
tableToPoints1Display.DataAxesGrid.YTitleFontFile = ''
tableToPoints1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
tableToPoints1Display.DataAxesGrid.ZTitleFontFile = ''
tableToPoints1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
tableToPoints1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
tableToPoints1Display.DataAxesGrid.XLabelFontFile = ''
tableToPoints1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
tableToPoints1Display.DataAxesGrid.YLabelFontFile = ''
tableToPoints1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
tableToPoints1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
tableToPoints1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
tableToPoints1Display.PolarAxes.PolarAxisTitleFontFile = ''
tableToPoints1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
tableToPoints1Display.PolarAxes.PolarAxisLabelFontFile = ''
tableToPoints1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
tableToPoints1Display.PolarAxes.LastRadialAxisTextFontFile = ''
tableToPoints1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
tableToPoints1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for kLUT in view renderView1
kLUTColorBar = GetScalarBar(kLUT, renderView1)
kLUTColorBar.Title = 'k'
kLUTColorBar.ComponentTitle = ''
kLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
kLUTColorBar.TitleFontFile = ''
kLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
kLUTColorBar.LabelFontFile = ''

# set color bar visibility
kLUTColorBar.Visibility = 1

# show color legend
tableToPoints1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'k'
kPWF = GetOpacityTransferFunction('k')
kPWF.Points = [0.0, 1.0, 0.5, 0.0, 6.793761253356934, 1.0, 0.5, 0.0, 8.928942680358887, 0.0, 0.5, 0.0, 17.760831832885742, 0.0, 0.5, 0.0, 19.701906204223633, 0.9924242496490479, 0.5, 0.0, 28.0, 1.0, 0.5, 0.0]
kPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(acsv)
# ----------------------------------------------------------------