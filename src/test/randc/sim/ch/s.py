# state file generated using paraview version 5.5.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.1

#### import the simple module from the paraview
import sys
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1998, 923]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5099512603467459, 0.49048581092909566, 0.5599339630758458]
renderView1.StereoType = 0
if 0:
    renderView1.CameraPosition = [2.813868569960714, -0.5074590100179297, 1.4878129079851992]
    renderView1.CameraFocalPoint = [0.5099512603467459, 0.49048581092909566, 0.5599339630758458]
    renderView1.CameraViewUp = [-0.31739082620003495, 0.13938956480260697, 0.9379944630264078]
    renderView1.CameraParallelScale = 0.6927889287942449
else:
    renderView1.CameraPosition = [0.39800458083951107, 0.6522102516464146, 1.0596291883480715]
    renderView1.CameraFocalPoint = [0.5100810080766677, 0.4898284971714025, 0.56021748483181]
    renderView1.CameraViewUp = [0.05872615649612774, 0.9531525463623923, -0.29673466583141267]
    renderView1.CameraParallelScale = 0.13897908169648868
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableOSPRay = 1
renderView1.Shadows = 1
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 5
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
partitfcsv = CSVReader(FileName=[sys.argv[2]])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=partitfcsv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Legacy VTK Reader'
s_vtk = LegacyVTKReader(FileNames=[sys.argv[1]])

# create a new 'Calculator'
calculator1 = Calculator(Input=tableToPoints1)
calculator1.ResultArrayName = 'cc'
calculator1.AttributeType = 'Point Data'
calculator1.Function = 'sin(1234567*c)'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from s_vtk
s_vtkDisplay = Show(s_vtk, renderView1)

# trace defaults for the display properties.
s_vtkDisplay.Representation = 'Surface'
s_vtkDisplay.AmbientColor = [0., 0., 0.]
s_vtkDisplay.ColorArrayName = [None, '']
s_vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s_vtkDisplay.SelectOrientationVectors = 'None'
s_vtkDisplay.ScaleFactor = 0.0801363028585911
s_vtkDisplay.SelectScaleArray = 'None'
s_vtkDisplay.GlyphType = 'Arrow'
s_vtkDisplay.GlyphTableIndexArray = 'None'
s_vtkDisplay.GaussianRadius = 0.00400681514292955
s_vtkDisplay.SetScaleArray = [None, '']
s_vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
s_vtkDisplay.OpacityArray = [None, '']
s_vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
s_vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
s_vtkDisplay.SelectionCellLabelFontFile = ''
s_vtkDisplay.SelectionPointLabelFontFile = ''
s_vtkDisplay.PolarAxes = 'PolarAxesRepresentation'
# change solid color
s_vtkDisplay.DiffuseColor = [0.6, 0.6, 0.6]

# change representation type
s_vtkDisplay.SetRepresentationType('Surface With Edges')

# Properties modified on s_0000vtkDisplay
s_vtkDisplay.RenderLinesAsTubes = 1

# Properties modified on s_0000vtkDisplay
s_vtkDisplay.LineWidth = 2.0

# Properties modified on s_0000vtkDisplay
s_vtkDisplay.EdgeColor = [0.3, 0.3, 0.3]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_vtkDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
s_vtkDisplay.DataAxesGrid.XTitleFontFile = ''
s_vtkDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
s_vtkDisplay.DataAxesGrid.YTitleFontFile = ''
s_vtkDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
s_vtkDisplay.DataAxesGrid.ZTitleFontFile = ''
s_vtkDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
s_vtkDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
s_vtkDisplay.DataAxesGrid.XLabelFontFile = ''
s_vtkDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
s_vtkDisplay.DataAxesGrid.YLabelFontFile = ''
s_vtkDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
s_vtkDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_vtkDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
s_vtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
s_vtkDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
s_vtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
s_vtkDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
s_vtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
s_vtkDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
s_vtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from calculator1
calculator1Display = Show(calculator1, renderView1)

# get color transfer function/color map for 'cc'
ccLUT = GetColorTransferFunction('cc')
ccLUT.RGBPoints = [-0.9999950949114859, 0.0, 0.0, 1.0, -0.6679959104213643, 0.0, 0.0, 1.0, -0.6659959153340744, 1.0, 0.0, 1.0, -0.33599672593124275, 1.0, 0.0, 1.0, -0.33399673084395265, 0.0, 1.0, 1.0, 2.4487334587819376e-06, 0.0, 1.0, 1.0, 0.002002443820748767, 0.0, 1.0, 0.0, 0.33200163322358056, 0.0, 1.0, 0.0, 0.3340016283108701, 1.0, 1.0, 0.0, 0.6640008177137019, 1.0, 1.0, 0.0, 0.6660008128009919, 1.0, 0.0, 0.0, 0.9999999923784034, 1.0, 0.0, 0.0]
ccLUT.ColorSpace = 'HSV'
ccLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.AmbientColor = [0.0, 0.0, 0.0]
calculator1Display.ColorArrayName = ['POINTS', 'cc']
calculator1Display.LookupTable = ccLUT
calculator1Display.PointSize = 5.0
calculator1Display.RenderPointsAsSpheres = 1
calculator1Display.OSPRayScaleArray = 'c'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 0.0804747700521705
calculator1Display.SelectScaleArray = 'None'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'None'
calculator1Display.GaussianRadius = 0.00402373850260852
calculator1Display.SetScaleArray = ['POINTS', 'c']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'c']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.SelectionCellLabelFontFile = ''
calculator1Display.SelectionPointLabelFontFile = ''
calculator1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [10027032.0, 0.0, 0.5, 0.0, 61035034.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [10027032.0, 0.0, 0.5, 0.0, 61035034.0, 1.0, 0.5, 0.0]

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

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'cc'
ccPWF = GetOpacityTransferFunction('cc')
ccPWF.Points = [-0.9999950949114859, 0.0, 0.5, 0.0, -0.7560932640841367, 0.0, 0.5, 0.0, -0.7560932640841367, 1.0, 0.5, 0.0, 0.64992909424709, 1.0, 0.5, 0.0, 0.64992909424709, 0.0, 0.5, 0.0, 0.9999999923784034, 0.0, 0.5, 0.0]
ccPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
#SetActiveSource(calculator1)
# ----------------------------------------------------------------

# reset view to fit data
renderView1.ResetCamera()


# save screenshot
SaveScreenshot('s.png', renderView1, ImageResolution=[1000, 1000])

