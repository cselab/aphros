# state file generated using paraview version 5.5.0-RC3

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.0-RC3

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Bar Chart View'
barChartView1 = CreateView('XYBarChartView')
barChartView1.ViewSize = [780, 861]
barChartView1.ChartTitleFontFile = ''
barChartView1.LegendPosition = [665, 819]
barChartView1.LeftAxisTitleFontFile = ''
barChartView1.LeftAxisRangeMaximum = 36.0
barChartView1.LeftAxisLabelFontFile = ''
barChartView1.BottomAxisTitleFontFile = ''
barChartView1.BottomAxisRangeMinimum = -1.0
barChartView1.BottomAxisRangeMaximum = 3.0
barChartView1.BottomAxisLabelFontFile = ''
barChartView1.RightAxisRangeMaximum = 6.66
barChartView1.RightAxisLabelFontFile = ''
barChartView1.TopAxisTitleFontFile = ''
barChartView1.TopAxisRangeMaximum = 6.66
barChartView1.TopAxisLabelFontFile = ''

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [780, 861]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.33833900510406495, 0.5902230089359284, 0.46418050678443906]
renderView1.StereoType = 0
renderView1.CameraPosition = [1.2168242160365017, 0.5193941462688718, 0.918918587775088]
renderView1.CameraFocalPoint = [0.7582231891897862, 0.5563693721978539, 0.6815288761517535]
renderView1.CameraViewUp = [0.29513792557831425, 0.8492019533454609, -0.43789227821415516]
renderView1.CameraParallelScale = 0.13480832062677817
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

# create a new 'Threshold'
threshold1 = Threshold(Input=None)
threshold1.Scalars = ['CELLS', 'u']
threshold1.ThresholdRange = [0.01, 0.99]

# create a new 'Legacy VTK Reader'
svtk = LegacyVTKReader(FileNames=['/home/kpetr/s/ch/src/test/advection/s.5.vtk'])

# create a new 'XML Structured Grid Reader'
pvts = XMLStructuredGridReader(FileName=['/home/kpetr/s/ch/src/test/advection/p.5.vts'])
pvts.CellArrayStatus = ['u', 'k', 'a', 'nx', 'ny', 'nz', 'kh', 'kp']

# create a new 'CSV Reader'
partit0_csv = CSVReader(FileName=['/home/kpetr/s/ch/src/test/advection/partit.0_0.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_1.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_2.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_3.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_4.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_5.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_6.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_7.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_8.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_9.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_10.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_11.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_12.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_13.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_14.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_15.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_16.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_17.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_18.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_19.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_20.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_21.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_22.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_23.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_24.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_25.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_26.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_27.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_28.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_29.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_30.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_31.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_32.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_33.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_34.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_35.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_36.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_37.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_38.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_39.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_40.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_41.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_42.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_43.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_44.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_45.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_46.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_47.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_48.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_49.csv', '/home/kpetr/s/ch/src/test/advection/partit.0_50.csv'])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=partit0_csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Calculator'
calculator1 = Calculator(Input=tableToPoints1)
calculator1.ResultArrayName = 'cc'
calculator1.Function = 'sin(c)'

# create a new 'Calculator'
calculator2 = Calculator(Input=threshold1)
calculator2.AttributeType = 'Cell Data'
calculator2.ResultArrayName = 'ke'
calculator2.Function = '2 / 0.0625'

# create a new 'Calculator'
calculator3 = Calculator(Input=calculator2)
calculator3.AttributeType = 'Cell Data'
calculator3.ResultArrayName = 'kpn'
calculator3.Function = 'kp/ke'

# create a new 'Calculator'
calculator4 = Calculator(Input=calculator3)
calculator4.AttributeType = 'Cell Data'
calculator4.ResultArrayName = 'khn'
calculator4.Function = 'kh/ke'

# create a new 'Histogram'
histogram1 = Histogram(Input=calculator4)
histogram1.SelectInputArray = ['CELLS', 'kpn']
histogram1.BinCount = 50
histogram1.UseCustomBinRanges = 1
histogram1.CustomBinRanges = [-1.0, 3.0]

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
histogram1Display.SeriesLabel = ['bin_extents', 'bin_extents', 'bin_values', 'curvature']
histogram1Display.SeriesColor = ['bin_extents', '0', '0', '0', 'bin_values', '0.889998', '0.100008', '0.110002']
histogram1Display.SeriesPlotCorner = ['bin_extents', '0', 'bin_values', '0']
histogram1Display.SeriesLabelPrefix = ''

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from svtk
svtkDisplay = Show(svtk, renderView1)

# trace defaults for the display properties.
svtkDisplay.Representation = 'Surface'
svtkDisplay.ColorArrayName = [None, '']
svtkDisplay.Opacity = 0.68
svtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
svtkDisplay.SelectOrientationVectors = 'None'
svtkDisplay.ScaleFactor = 0.0387654006481171
svtkDisplay.SelectScaleArray = 'None'
svtkDisplay.GlyphType = 'Arrow'
svtkDisplay.GlyphTableIndexArray = 'None'
svtkDisplay.GaussianRadius = 0.00193827003240585
svtkDisplay.SetScaleArray = [None, '']
svtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
svtkDisplay.OpacityArray = [None, '']
svtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
svtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
svtkDisplay.SelectionCellLabelFontFile = ''
svtkDisplay.SelectionPointLabelFontFile = ''
svtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
svtkDisplay.DataAxesGrid.XTitleFontFile = ''
svtkDisplay.DataAxesGrid.YTitleFontFile = ''
svtkDisplay.DataAxesGrid.ZTitleFontFile = ''
svtkDisplay.DataAxesGrid.XLabelFontFile = ''
svtkDisplay.DataAxesGrid.YLabelFontFile = ''
svtkDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
svtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
svtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
svtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
svtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from calculator1
calculator1Display = Show(calculator1, renderView1)

# get color transfer function/color map for 'cc'
ccLUT = GetColorTransferFunction('cc')
ccLUT.EnableOpacityMapping = 1
ccLUT.RGBPoints = [-0.999999046236321, 0.0, 0.0, 1.0, -0.66799941495941, 0.0, 0.0, 1.0, -0.665999417180633, 1.0, 0.0, 1.0, -0.335999783682499, 1.0, 0.0, 1.0, -0.333999785903722, 0.0, 1.0, 1.0, -1.56848034560397e-07, 0.0, 1.0, 1.0, 0.00199984093074201, 0.0, 1.0, 0.0, 0.331999474428876, 0.0, 1.0, 0.0, 0.333999472207653, 1.0, 1.0, 0.0, 0.663999105705787, 1.0, 1.0, 0.0, 0.665999103484564, 1.0, 0.0, 0.0, 0.999998732540252, 1.0, 0.0, 0.0]
ccLUT.ColorSpace = 'HSV'
ccLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
calculator1Display.Representation = 'Point Gaussian'
calculator1Display.ColorArrayName = ['POINTS', 'cc']
calculator1Display.LookupTable = ccLUT
calculator1Display.OSPRayScaleArray = 'cc'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 0.0394942
calculator1Display.SelectScaleArray = 'cc'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'cc'
calculator1Display.GaussianRadius = 0.0018562274
calculator1Display.SetScaleArray = ['POINTS', 'cc']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'cc']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.SelectionCellLabelFontFile = ''
calculator1Display.SelectionPointLabelFontFile = ''
calculator1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [-0.9992737805272319, 0.0, 0.5, 0.0, 0.9936602925070573, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [-0.9992737805272319, 0.0, 0.5, 0.0, 0.9936602925070573, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calculator1Display.DataAxesGrid.XTitleFontFile = ''
calculator1Display.DataAxesGrid.YTitleFontFile = ''
calculator1Display.DataAxesGrid.ZTitleFontFile = ''
calculator1Display.DataAxesGrid.XLabelFontFile = ''
calculator1Display.DataAxesGrid.YLabelFontFile = ''
calculator1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calculator1Display.PolarAxes.PolarAxisTitleFontFile = ''
calculator1Display.PolarAxes.PolarAxisLabelFontFile = ''
calculator1Display.PolarAxes.LastRadialAxisTextFontFile = ''
calculator1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for ccLUT in view renderView1
ccLUTColorBar = GetScalarBar(ccLUT, renderView1)
ccLUTColorBar.Title = 'cc'
ccLUTColorBar.ComponentTitle = ''
ccLUTColorBar.TitleFontFile = ''
ccLUTColorBar.LabelFontFile = ''

# set color bar visibility
ccLUTColorBar.Visibility = 1

# show color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'cc'
ccPWF = GetOpacityTransferFunction('cc')
ccPWF.Points = [-0.999999046236321, 0.0, 0.5, 0.0, 0.0848326911321762, 0.0, 0.5, 0.0, 0.0848326911321762, 1.0, 0.5, 0.0, 0.24935694346773, 0.993421077728271, 0.5, 0.0, 0.24935694346773, 0.0, 0.5, 0.0, 0.999998732540252, 0.0, 0.5, 0.0]
ccPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(calculator4)
# ----------------------------------------------------------------