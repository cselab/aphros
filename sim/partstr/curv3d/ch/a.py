# state file generated using paraview version 5.5.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1769, 1205]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5100770108401775, 0.47003214433789253, 0.559935986995697]
renderView1.StereoType = 0
renderView1.CameraPosition = [2.6175173876539772, 1.9956633671894988, 1.2075107981789321]
renderView1.CameraFocalPoint = [0.510077010840177, 0.47003214433789225, 0.5599359869956971]
renderView1.CameraViewUp = [0.1944234025884112, -0.5984618922823668, 0.7772019711836515]
renderView1.CameraParallelScale = 0.6939154699620841
renderView1.Background = [1.0, 1.0, 1.0]
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
partit0csv = CSVReader(FileName=['/home/kpetr/s/ch/sim/partstr/curv3d/ch/partit.0.csv'])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=partit0csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Calculator'
calculator1 = Calculator(Input=tableToPoints1)
calculator1.ResultArrayName = 'cc'
calculator1.Function = 'sin(1237*c)'

# create a new 'Legacy VTK Reader'
s0vtk = LegacyVTKReader(FileNames=['/home/kpetr/s/ch/sim/partstr/curv3d/ch/s.0.vtk'])

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from s0vtk
s0vtkDisplay = Show(s0vtk, renderView1)

# trace defaults for the display properties.
s0vtkDisplay.Representation = 'Surface'
s0vtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
s0vtkDisplay.ColorArrayName = [None, '']
s0vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s0vtkDisplay.SelectOrientationVectors = 'None'
s0vtkDisplay.ScaleFactor = 0.08013630285859108
s0vtkDisplay.SelectScaleArray = 'None'
s0vtkDisplay.GlyphType = 'Arrow'
s0vtkDisplay.GlyphTableIndexArray = 'None'
s0vtkDisplay.GaussianRadius = 0.004006815142929554
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

# show data from calculator1
calculator1Display = Show(calculator1, renderView1)

# get color transfer function/color map for 'cc'
ccLUT = GetColorTransferFunction('cc')
ccLUT.AutomaticRescaleRangeMode = 'Never'
ccLUT.EnableOpacityMapping = 1
ccLUT.RGBPoints = [0.0, 0.0, 0.0, 1.0, 0.0166, 0.0, 0.0, 1.0, 0.016700000000000003, 1.0, 0.0, 1.0, 0.0332, 1.0, 0.0, 1.0, 0.0333, 0.0, 1.0, 1.0, 0.05, 0.0, 1.0, 1.0, 0.050100000000000006, 0.0, 1.0, 0.0, 0.0666, 0.0, 1.0, 0.0, 0.06670000000000001, 1.0, 1.0, 0.0, 0.0832, 1.0, 1.0, 0.0, 0.0833, 1.0, 0.0, 0.0, 0.1, 1.0, 0.0, 0.0]
ccLUT.ColorSpace = 'HSV'
ccLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.AmbientColor = [0.0, 0.0, 0.0]
calculator1Display.ColorArrayName = ['POINTS', 'cc']
calculator1Display.LookupTable = ccLUT
calculator1Display.PointSize = 8.0
calculator1Display.RenderPointsAsSpheres = 1
calculator1Display.OSPRayScaleArray = 'c'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 0.08047477005217046
calculator1Display.SelectScaleArray = 'None'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'None'
calculator1Display.GaussianRadius = 0.004023738502608522
calculator1Display.SetScaleArray = ['POINTS', 'c']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'c']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.SelectionCellLabelFontFile = ''
calculator1Display.SelectionPointLabelFontFile = ''
calculator1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [5010013.0, 0.0, 0.5, 0.0, 30019017.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [5010013.0, 0.0, 0.5, 0.0, 30019017.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calculator1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
calculator1Display.DataAxesGrid.XTitleFontFile = ''
calculator1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
calculator1Display.DataAxesGrid.YTitleFontFile = ''
calculator1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
calculator1Display.DataAxesGrid.ZTitleFontFile = ''
calculator1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
calculator1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
calculator1Display.DataAxesGrid.XLabelFontFile = ''
calculator1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
calculator1Display.DataAxesGrid.YLabelFontFile = ''
calculator1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
calculator1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calculator1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
calculator1Display.PolarAxes.PolarAxisTitleFontFile = ''
calculator1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
calculator1Display.PolarAxes.PolarAxisLabelFontFile = ''
calculator1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
calculator1Display.PolarAxes.LastRadialAxisTextFontFile = ''
calculator1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
calculator1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for ccLUT in view renderView1
ccLUTColorBar = GetScalarBar(ccLUT, renderView1)
ccLUTColorBar.Title = 'cc'
ccLUTColorBar.ComponentTitle = ''
ccLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
ccLUTColorBar.TitleFontFile = ''
ccLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
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
ccPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.012195121496915817, 0.0, 0.5, 0.0, 0.012195121496915817, 1.0, 0.5, 0.0, 0.08249641209840775, 1.0, 0.5, 0.0, 0.08249641209840775, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0]
ccPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(partit0csv)
# ----------------------------------------------------------------