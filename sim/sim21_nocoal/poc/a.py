# state file generated using paraview version 5.6.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [583, 400]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 0.5, 0.015625]
renderView1.UseLight = 0
renderView1.StereoType = 0
renderView1.CameraPosition = [0.5, 0.5, 2.7483427307589823]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.015625]
renderView1.CameraParallelScale = 0.565719768687297
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelFontFile = ''

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [583, 400]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [0.5, 0.5, 0.015625]
renderView2.UseLight = 0
renderView2.StereoType = 0
renderView2.CameraPosition = [0.5, 0.5, 2.7483427307589823]
renderView2.CameraFocalPoint = [0.5, 0.5, 0.015625]
renderView2.CameraParallelScale = 0.565719768687297
renderView2.CameraParallelProjection = 1
renderView2.Background = [0.0, 0.0, 0.0]
renderView2.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView2.AxesGrid.XTitleFontFile = ''
renderView2.AxesGrid.YTitleFontFile = ''
renderView2.AxesGrid.ZTitleFontFile = ''
renderView2.AxesGrid.XLabelFontFile = ''
renderView2.AxesGrid.YLabelFontFile = ''
renderView2.AxesGrid.ZLabelFontFile = ''

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.ViewSize = [583, 400]
renderView3.AxesGrid = 'GridAxes3DActor'
renderView3.OrientationAxesVisibility = 0
renderView3.CenterOfRotation = [0.5, 0.5, 0.015625]
renderView3.UseLight = 0
renderView3.StereoType = 0
renderView3.CameraPosition = [0.5, 0.5, 2.7483427307589823]
renderView3.CameraFocalPoint = [0.5, 0.5, 0.015625]
renderView3.CameraParallelScale = 0.565719768687297
renderView3.CameraParallelProjection = 1
renderView3.Background = [0.0, 0.0, 0.0]
renderView3.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView3.AxesGrid.XTitleFontFile = ''
renderView3.AxesGrid.YTitleFontFile = ''
renderView3.AxesGrid.ZTitleFontFile = ''
renderView3.AxesGrid.XLabelFontFile = ''
renderView3.AxesGrid.YLabelFontFile = ''
renderView3.AxesGrid.ZLabelFontFile = ''

# Create a new 'Render View'
renderView4 = CreateView('RenderView')
renderView4.ViewSize = [583, 400]
renderView4.AxesGrid = 'GridAxes3DActor'
renderView4.OrientationAxesVisibility = 0
renderView4.CenterOfRotation = [0.5, 0.5, 0.015625]
renderView4.UseLight = 0
renderView4.StereoType = 0
renderView4.CameraPosition = [0.5, 0.5, 2.7483427307589823]
renderView4.CameraFocalPoint = [0.5, 0.5, 0.015625]
renderView4.CameraParallelScale = 0.565719768687297
renderView4.CameraParallelProjection = 1
renderView4.Background = [0.0, 0.0, 0.0]
renderView4.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView4.AxesGrid.XTitleFontFile = ''
renderView4.AxesGrid.YTitleFontFile = ''
renderView4.AxesGrid.ZTitleFontFile = ''
renderView4.AxesGrid.XLabelFontFile = ''
renderView4.AxesGrid.YLabelFontFile = ''
renderView4.AxesGrid.ZLabelFontFile = ''

# Create a new 'Render View'
renderView5 = CreateView('RenderView')
renderView5.ViewSize = [1175, 830]
renderView5.InteractionMode = '2D'
renderView5.AxesGrid = 'GridAxes3DActor'
renderView5.OrientationAxesVisibility = 0
renderView5.CenterOfRotation = [0.5098712518811226, 0.4999999850988388, 0.015625]
renderView5.StereoType = 0
renderView5.CameraPosition = [0.5098712518811226, 0.4999999850988388, 1.7828294314862205]
renderView5.CameraFocalPoint = [0.5098712518811226, 0.4999999850988388, 0.015625]
renderView5.CameraParallelScale = 0.5534372577844301
renderView5.Background = [0.0, 0.0, 0.0]
renderView5.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView5.AxesGrid.XTitleFontFile = ''
renderView5.AxesGrid.YTitleFontFile = ''
renderView5.AxesGrid.ZTitleFontFile = ''
renderView5.AxesGrid.XLabelFontFile = ''
renderView5.AxesGrid.YLabelFontFile = ''
renderView5.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView5)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
s_vtk = LegacyVTKReader(FileNames=['/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0000.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0001.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0002.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0003.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0004.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0005.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0006.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0007.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0008.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0009.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0010.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0011.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0012.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0013.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0014.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0015.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0016.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0017.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0018.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0019.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0020.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0021.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0022.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0023.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0024.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0025.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0026.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0027.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0028.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0029.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0030.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0031.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0032.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0033.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0034.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0035.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0036.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0037.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0038.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0039.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0040.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0041.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0042.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0043.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0044.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0045.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0046.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0047.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0048.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0049.vtk', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/s_0050.vtk'])

# create a new 'XML Structured Grid Reader'
p_vts = XMLStructuredGridReader(FileName=['/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0000.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0001.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0002.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0003.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0004.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0005.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0006.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0007.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0008.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0009.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0010.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0011.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0012.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0013.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0014.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0015.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0016.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0017.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0018.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0019.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0020.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0021.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0022.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0023.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0024.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0025.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0026.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0027.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0028.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0029.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0030.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0031.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0032.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0033.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0034.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0035.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0036.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0037.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0038.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0039.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0040.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0041.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0042.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0043.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0044.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0045.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0046.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0047.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0048.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0049.vts', '/home/petr/work/ethz/proj/electrochem/chtmp/drain/Jul24/p_0050.vts'])
p_vts.CellArrayStatus = ['p', 'vf', 'sig', 'vx', 'vy', 'nx', 'ny', 'n1x', 'n1y', 'n2x', 'n2y', 'vf1', 'vf2', 'cl', 'mask1', 'mask2', 'mul', 'mul2', 'cl1', 'cl2', 'dep1', 'dep2']

# create a new 'Outline'
outline1 = Outline(Input=p_vts)

# create a new 'Threshold'
threshold1 = Threshold(Input=p_vts)
threshold1.Scalars = ['CELLS', 'cl0']
threshold1.ThresholdRange = [0.0, 11.0]

# create a new 'Extract Edges'
extractEdges1 = ExtractEdges(Input=p_vts)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from p_vts
p_vtsDisplay = Show(p_vts, renderView1)

# get color transfer function/color map for 'vf1'
vf1LUT = GetColorTransferFunction('vf1')
vf1LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.5000000000000001, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]
vf1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'vf1'
vf1PWF = GetOpacityTransferFunction('vf1')
vf1PWF.Points = [0.0, 1.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
vf1PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
p_vtsDisplay.Representation = 'Surface With Edges'
p_vtsDisplay.ColorArrayName = ['CELLS', 'vf1']
p_vtsDisplay.LookupTable = vf1LUT
p_vtsDisplay.Opacity = 0.5
p_vtsDisplay.EdgeColor = [0.3137254901960784, 0.3137254901960784, 0.3137254901960784]
p_vtsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
p_vtsDisplay.SelectOrientationVectors = 'None'
p_vtsDisplay.ScaleFactor = 0.1
p_vtsDisplay.SelectScaleArray = 'None'
p_vtsDisplay.GlyphType = 'Arrow'
p_vtsDisplay.GlyphTableIndexArray = 'None'
p_vtsDisplay.GaussianRadius = 0.005
p_vtsDisplay.SetScaleArray = [None, '']
p_vtsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
p_vtsDisplay.OpacityArray = [None, '']
p_vtsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
p_vtsDisplay.DataAxesGrid = 'GridAxesRepresentation'
p_vtsDisplay.SelectionCellLabelFontFile = ''
p_vtsDisplay.SelectionPointLabelFontFile = ''
p_vtsDisplay.PolarAxes = 'PolarAxesRepresentation'
p_vtsDisplay.ScalarOpacityFunction = vf1PWF
p_vtsDisplay.ScalarOpacityUnitDistance = 0.08839374228030177

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
p_vtsDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
p_vtsDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
p_vtsDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
p_vtsDisplay.DataAxesGrid.XTitleFontFile = ''
p_vtsDisplay.DataAxesGrid.YTitleFontFile = ''
p_vtsDisplay.DataAxesGrid.ZTitleFontFile = ''
p_vtsDisplay.DataAxesGrid.XLabelFontFile = ''
p_vtsDisplay.DataAxesGrid.YLabelFontFile = ''
p_vtsDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
p_vtsDisplay.PolarAxes.PolarAxisTitleFontFile = ''
p_vtsDisplay.PolarAxes.PolarAxisLabelFontFile = ''
p_vtsDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
p_vtsDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from s_vtk
s_vtkDisplay = Show(s_vtk, renderView1)

# get color transfer function/color map for 'l'
lLUT = GetColorTransferFunction('l')
lLUT.RGBPoints = [0.0, 0.0, 0.0, 1.0, 0.664, 0.0, 0.0, 1.0, 0.668, 1.0, 0.0, 1.0, 1.328, 1.0, 0.0, 1.0, 1.332, 0.0, 1.0, 1.0, 2.0, 0.0, 1.0, 1.0, 2.004, 0.0, 1.0, 0.0, 2.664, 0.0, 1.0, 0.0, 2.6680000000000006, 1.0, 1.0, 0.0, 3.328, 1.0, 1.0, 0.0, 3.3319999999999994, 1.0, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0]
lLUT.ColorSpace = 'HSV'
lLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
s_vtkDisplay.Representation = 'Wireframe'
s_vtkDisplay.ColorArrayName = ['CELLS', 'l']
s_vtkDisplay.LookupTable = lLUT
s_vtkDisplay.PointSize = 30.0
s_vtkDisplay.LineWidth = 4.0
s_vtkDisplay.RenderPointsAsSpheres = 1
s_vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s_vtkDisplay.SelectOrientationVectors = 'None'
s_vtkDisplay.ScaleFactor = 0.02587677538394928
s_vtkDisplay.SelectScaleArray = 'c'
s_vtkDisplay.GlyphType = 'Arrow'
s_vtkDisplay.GlyphTableIndexArray = 'c'
s_vtkDisplay.GaussianRadius = 0.001293838769197464
s_vtkDisplay.SetScaleArray = [None, '']
s_vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
s_vtkDisplay.OpacityArray = [None, '']
s_vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
s_vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
s_vtkDisplay.SelectionCellLabelFontFile = ''
s_vtkDisplay.SelectionPointLabelFontFile = ''
s_vtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_vtkDisplay.DataAxesGrid.XTitleFontFile = ''
s_vtkDisplay.DataAxesGrid.YTitleFontFile = ''
s_vtkDisplay.DataAxesGrid.ZTitleFontFile = ''
s_vtkDisplay.DataAxesGrid.XLabelFontFile = ''
s_vtkDisplay.DataAxesGrid.YLabelFontFile = ''
s_vtkDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_vtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
s_vtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
s_vtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
s_vtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from outline1
outline1Display = Show(outline1, renderView1)

# trace defaults for the display properties.
outline1Display.Representation = 'Surface'
outline1Display.ColorArrayName = [None, '']
outline1Display.OSPRayScaleFunction = 'PiecewiseFunction'
outline1Display.SelectOrientationVectors = 'None'
outline1Display.ScaleFactor = 0.1
outline1Display.SelectScaleArray = 'None'
outline1Display.GlyphType = 'Arrow'
outline1Display.GlyphTableIndexArray = 'None'
outline1Display.GaussianRadius = 0.005
outline1Display.SetScaleArray = [None, '']
outline1Display.ScaleTransferFunction = 'PiecewiseFunction'
outline1Display.OpacityArray = [None, '']
outline1Display.OpacityTransferFunction = 'PiecewiseFunction'
outline1Display.DataAxesGrid = 'GridAxesRepresentation'
outline1Display.SelectionCellLabelFontFile = ''
outline1Display.SelectionPointLabelFontFile = ''
outline1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
outline1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
outline1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
outline1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
outline1Display.DataAxesGrid.XTitleFontFile = ''
outline1Display.DataAxesGrid.YTitleFontFile = ''
outline1Display.DataAxesGrid.ZTitleFontFile = ''
outline1Display.DataAxesGrid.XLabelFontFile = ''
outline1Display.DataAxesGrid.YLabelFontFile = ''
outline1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
outline1Display.PolarAxes.PolarAxisTitleFontFile = ''
outline1Display.PolarAxes.PolarAxisLabelFontFile = ''
outline1Display.PolarAxes.LastRadialAxisTextFontFile = ''
outline1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for vf1LUT in view renderView1
vf1LUTColorBar = GetScalarBar(vf1LUT, renderView1)
vf1LUTColorBar.Position = [0.8028335301062575, 0.07478260869565218]
vf1LUTColorBar.Title = 'vf1'
vf1LUTColorBar.ComponentTitle = ''
vf1LUTColorBar.HorizontalTitle = 1
vf1LUTColorBar.TitleFontFamily = 'Courier'
vf1LUTColorBar.TitleFontFile = ''
vf1LUTColorBar.TitleBold = 1
vf1LUTColorBar.TitleFontSize = 20
vf1LUTColorBar.LabelFontFile = ''
vf1LUTColorBar.AutomaticLabelFormat = 0
vf1LUTColorBar.LabelFormat = '%-#6.1f'
vf1LUTColorBar.AddRangeLabels = 0
vf1LUTColorBar.RangeLabelFormat = '%-#6.1f'

# set color bar visibility
vf1LUTColorBar.Visibility = 1

# show color legend
p_vtsDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from s_vtk
s_vtkDisplay_1 = Show(s_vtk, renderView2)

# trace defaults for the display properties.
s_vtkDisplay_1.Representation = 'Wireframe'
s_vtkDisplay_1.ColorArrayName = ['CELLS', 'l']
s_vtkDisplay_1.LookupTable = lLUT
s_vtkDisplay_1.Opacity = 0.75
s_vtkDisplay_1.PointSize = 30.0
s_vtkDisplay_1.LineWidth = 4.0
s_vtkDisplay_1.RenderPointsAsSpheres = 1
s_vtkDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
s_vtkDisplay_1.SelectOrientationVectors = 'None'
s_vtkDisplay_1.ScaleFactor = 0.02587677538394928
s_vtkDisplay_1.SelectScaleArray = 'c'
s_vtkDisplay_1.GlyphType = 'Arrow'
s_vtkDisplay_1.GlyphTableIndexArray = 'c'
s_vtkDisplay_1.GaussianRadius = 0.001293838769197464
s_vtkDisplay_1.SetScaleArray = [None, '']
s_vtkDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
s_vtkDisplay_1.OpacityArray = [None, '']
s_vtkDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
s_vtkDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
s_vtkDisplay_1.SelectionCellLabelFontFile = ''
s_vtkDisplay_1.SelectionPointLabelFontFile = ''
s_vtkDisplay_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_vtkDisplay_1.DataAxesGrid.XTitleFontFile = ''
s_vtkDisplay_1.DataAxesGrid.YTitleFontFile = ''
s_vtkDisplay_1.DataAxesGrid.ZTitleFontFile = ''
s_vtkDisplay_1.DataAxesGrid.XLabelFontFile = ''
s_vtkDisplay_1.DataAxesGrid.YLabelFontFile = ''
s_vtkDisplay_1.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_vtkDisplay_1.PolarAxes.PolarAxisTitleFontFile = ''
s_vtkDisplay_1.PolarAxes.PolarAxisLabelFontFile = ''
s_vtkDisplay_1.PolarAxes.LastRadialAxisTextFontFile = ''
s_vtkDisplay_1.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from p_vts
p_vtsDisplay_1 = Show(p_vts, renderView2)

# get color transfer function/color map for 'vf2'
vf2LUT = GetColorTransferFunction('vf2')
vf2LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.0, 0.705882, 0.0156863, 0.14902]
vf2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'vf2'
vf2PWF = GetOpacityTransferFunction('vf2')
vf2PWF.Points = [0.0, 1.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
vf2PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
p_vtsDisplay_1.Representation = 'Surface With Edges'
p_vtsDisplay_1.ColorArrayName = ['CELLS', 'vf2']
p_vtsDisplay_1.LookupTable = vf2LUT
p_vtsDisplay_1.Opacity = 0.5
p_vtsDisplay_1.EdgeColor = [0.3137254901960784, 0.3137254901960784, 0.3137254901960784]
p_vtsDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
p_vtsDisplay_1.SelectOrientationVectors = 'None'
p_vtsDisplay_1.ScaleFactor = 0.1
p_vtsDisplay_1.SelectScaleArray = 'None'
p_vtsDisplay_1.GlyphType = 'Arrow'
p_vtsDisplay_1.GlyphTableIndexArray = 'None'
p_vtsDisplay_1.GaussianRadius = 0.005
p_vtsDisplay_1.SetScaleArray = [None, '']
p_vtsDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
p_vtsDisplay_1.OpacityArray = [None, '']
p_vtsDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
p_vtsDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
p_vtsDisplay_1.SelectionCellLabelFontFile = ''
p_vtsDisplay_1.SelectionPointLabelFontFile = ''
p_vtsDisplay_1.PolarAxes = 'PolarAxesRepresentation'
p_vtsDisplay_1.ScalarOpacityFunction = vf2PWF
p_vtsDisplay_1.ScalarOpacityUnitDistance = 0.08839374228030177

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
p_vtsDisplay_1.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
p_vtsDisplay_1.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
p_vtsDisplay_1.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
p_vtsDisplay_1.DataAxesGrid.XTitleFontFile = ''
p_vtsDisplay_1.DataAxesGrid.YTitleFontFile = ''
p_vtsDisplay_1.DataAxesGrid.ZTitleFontFile = ''
p_vtsDisplay_1.DataAxesGrid.XLabelFontFile = ''
p_vtsDisplay_1.DataAxesGrid.YLabelFontFile = ''
p_vtsDisplay_1.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
p_vtsDisplay_1.PolarAxes.PolarAxisTitleFontFile = ''
p_vtsDisplay_1.PolarAxes.PolarAxisLabelFontFile = ''
p_vtsDisplay_1.PolarAxes.LastRadialAxisTextFontFile = ''
p_vtsDisplay_1.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from outline1
outline1Display_1 = Show(outline1, renderView2)

# trace defaults for the display properties.
outline1Display_1.Representation = 'Surface'
outline1Display_1.ColorArrayName = [None, '']
outline1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
outline1Display_1.SelectOrientationVectors = 'None'
outline1Display_1.ScaleFactor = 0.1
outline1Display_1.SelectScaleArray = 'None'
outline1Display_1.GlyphType = 'Arrow'
outline1Display_1.GlyphTableIndexArray = 'None'
outline1Display_1.GaussianRadius = 0.005
outline1Display_1.SetScaleArray = [None, '']
outline1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
outline1Display_1.OpacityArray = [None, '']
outline1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
outline1Display_1.DataAxesGrid = 'GridAxesRepresentation'
outline1Display_1.SelectionCellLabelFontFile = ''
outline1Display_1.SelectionPointLabelFontFile = ''
outline1Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
outline1Display_1.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
outline1Display_1.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
outline1Display_1.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
outline1Display_1.DataAxesGrid.XTitleFontFile = ''
outline1Display_1.DataAxesGrid.YTitleFontFile = ''
outline1Display_1.DataAxesGrid.ZTitleFontFile = ''
outline1Display_1.DataAxesGrid.XLabelFontFile = ''
outline1Display_1.DataAxesGrid.YLabelFontFile = ''
outline1Display_1.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
outline1Display_1.PolarAxes.PolarAxisTitleFontFile = ''
outline1Display_1.PolarAxes.PolarAxisLabelFontFile = ''
outline1Display_1.PolarAxes.LastRadialAxisTextFontFile = ''
outline1Display_1.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for vf2LUT in view renderView2
vf2LUTColorBar = GetScalarBar(vf2LUT, renderView2)
vf2LUTColorBar.Position = [0.8087367178276268, 0.07304347826086952]
vf2LUTColorBar.Title = 'vf2'
vf2LUTColorBar.ComponentTitle = ''
vf2LUTColorBar.HorizontalTitle = 1
vf2LUTColorBar.TitleFontFamily = 'Courier'
vf2LUTColorBar.TitleFontFile = ''
vf2LUTColorBar.TitleBold = 1
vf2LUTColorBar.TitleFontSize = 20
vf2LUTColorBar.LabelFontFile = ''
vf2LUTColorBar.AutomaticLabelFormat = 0
vf2LUTColorBar.LabelFormat = '%-#6.1f'
vf2LUTColorBar.AddRangeLabels = 0
vf2LUTColorBar.RangeLabelFormat = '%-#6.1f'
vf2LUTColorBar.ScalarBarLength = 0.32999999999999996

# set color bar visibility
vf2LUTColorBar.Visibility = 1

# show color legend
p_vtsDisplay_1.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView3'
# ----------------------------------------------------------------

# show data from p_vts
p_vtsDisplay_2 = Show(p_vts, renderView3)

# get color transfer function/color map for 'mask2'
mask2LUT = GetColorTransferFunction('mask2')
mask2LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.5000000000000001, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]
mask2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'mask2'
mask2PWF = GetOpacityTransferFunction('mask2')
mask2PWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
mask2PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
p_vtsDisplay_2.Representation = 'Surface With Edges'
p_vtsDisplay_2.ColorArrayName = ['CELLS', 'mask2']
p_vtsDisplay_2.LookupTable = mask2LUT
p_vtsDisplay_2.Opacity = 0.5
p_vtsDisplay_2.EdgeColor = [0.3137254901960784, 0.3137254901960784, 0.3137254901960784]
p_vtsDisplay_2.OSPRayScaleFunction = 'PiecewiseFunction'
p_vtsDisplay_2.SelectOrientationVectors = 'None'
p_vtsDisplay_2.ScaleFactor = 0.1
p_vtsDisplay_2.SelectScaleArray = 'None'
p_vtsDisplay_2.GlyphType = 'Arrow'
p_vtsDisplay_2.GlyphTableIndexArray = 'None'
p_vtsDisplay_2.GaussianRadius = 0.005
p_vtsDisplay_2.SetScaleArray = [None, '']
p_vtsDisplay_2.ScaleTransferFunction = 'PiecewiseFunction'
p_vtsDisplay_2.OpacityArray = [None, '']
p_vtsDisplay_2.OpacityTransferFunction = 'PiecewiseFunction'
p_vtsDisplay_2.DataAxesGrid = 'GridAxesRepresentation'
p_vtsDisplay_2.SelectionCellLabelFontFile = ''
p_vtsDisplay_2.SelectionPointLabelFontFile = ''
p_vtsDisplay_2.PolarAxes = 'PolarAxesRepresentation'
p_vtsDisplay_2.ScalarOpacityFunction = mask2PWF
p_vtsDisplay_2.ScalarOpacityUnitDistance = 0.08839374228030177

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
p_vtsDisplay_2.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
p_vtsDisplay_2.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
p_vtsDisplay_2.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
p_vtsDisplay_2.DataAxesGrid.XTitleFontFile = ''
p_vtsDisplay_2.DataAxesGrid.YTitleFontFile = ''
p_vtsDisplay_2.DataAxesGrid.ZTitleFontFile = ''
p_vtsDisplay_2.DataAxesGrid.XLabelFontFile = ''
p_vtsDisplay_2.DataAxesGrid.YLabelFontFile = ''
p_vtsDisplay_2.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
p_vtsDisplay_2.PolarAxes.PolarAxisTitleFontFile = ''
p_vtsDisplay_2.PolarAxes.PolarAxisLabelFontFile = ''
p_vtsDisplay_2.PolarAxes.LastRadialAxisTextFontFile = ''
p_vtsDisplay_2.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from s_vtk
s_vtkDisplay_2 = Show(s_vtk, renderView3)

# trace defaults for the display properties.
s_vtkDisplay_2.Representation = 'Wireframe'
s_vtkDisplay_2.ColorArrayName = ['CELLS', 'l']
s_vtkDisplay_2.LookupTable = lLUT
s_vtkDisplay_2.PointSize = 30.0
s_vtkDisplay_2.LineWidth = 4.0
s_vtkDisplay_2.RenderPointsAsSpheres = 1
s_vtkDisplay_2.OSPRayScaleFunction = 'PiecewiseFunction'
s_vtkDisplay_2.SelectOrientationVectors = 'None'
s_vtkDisplay_2.ScaleFactor = 0.02587677538394928
s_vtkDisplay_2.SelectScaleArray = 'c'
s_vtkDisplay_2.GlyphType = 'Arrow'
s_vtkDisplay_2.GlyphTableIndexArray = 'c'
s_vtkDisplay_2.GaussianRadius = 0.001293838769197464
s_vtkDisplay_2.SetScaleArray = [None, '']
s_vtkDisplay_2.ScaleTransferFunction = 'PiecewiseFunction'
s_vtkDisplay_2.OpacityArray = [None, '']
s_vtkDisplay_2.OpacityTransferFunction = 'PiecewiseFunction'
s_vtkDisplay_2.DataAxesGrid = 'GridAxesRepresentation'
s_vtkDisplay_2.SelectionCellLabelFontFile = ''
s_vtkDisplay_2.SelectionPointLabelFontFile = ''
s_vtkDisplay_2.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_vtkDisplay_2.DataAxesGrid.XTitleFontFile = ''
s_vtkDisplay_2.DataAxesGrid.YTitleFontFile = ''
s_vtkDisplay_2.DataAxesGrid.ZTitleFontFile = ''
s_vtkDisplay_2.DataAxesGrid.XLabelFontFile = ''
s_vtkDisplay_2.DataAxesGrid.YLabelFontFile = ''
s_vtkDisplay_2.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_vtkDisplay_2.PolarAxes.PolarAxisTitleFontFile = ''
s_vtkDisplay_2.PolarAxes.PolarAxisLabelFontFile = ''
s_vtkDisplay_2.PolarAxes.LastRadialAxisTextFontFile = ''
s_vtkDisplay_2.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from outline1
outline1Display_2 = Show(outline1, renderView3)

# trace defaults for the display properties.
outline1Display_2.Representation = 'Wireframe'
outline1Display_2.ColorArrayName = [None, '']
outline1Display_2.OSPRayScaleFunction = 'PiecewiseFunction'
outline1Display_2.SelectOrientationVectors = 'None'
outline1Display_2.ScaleFactor = 0.1
outline1Display_2.SelectScaleArray = 'None'
outline1Display_2.GlyphType = 'Arrow'
outline1Display_2.GlyphTableIndexArray = 'None'
outline1Display_2.GaussianRadius = 0.005
outline1Display_2.SetScaleArray = [None, '']
outline1Display_2.ScaleTransferFunction = 'PiecewiseFunction'
outline1Display_2.OpacityArray = [None, '']
outline1Display_2.OpacityTransferFunction = 'PiecewiseFunction'
outline1Display_2.DataAxesGrid = 'GridAxesRepresentation'
outline1Display_2.SelectionCellLabelFontFile = ''
outline1Display_2.SelectionPointLabelFontFile = ''
outline1Display_2.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
outline1Display_2.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
outline1Display_2.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
outline1Display_2.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
outline1Display_2.DataAxesGrid.XTitleFontFile = ''
outline1Display_2.DataAxesGrid.YTitleFontFile = ''
outline1Display_2.DataAxesGrid.ZTitleFontFile = ''
outline1Display_2.DataAxesGrid.XLabelFontFile = ''
outline1Display_2.DataAxesGrid.YLabelFontFile = ''
outline1Display_2.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
outline1Display_2.PolarAxes.PolarAxisTitleFontFile = ''
outline1Display_2.PolarAxes.PolarAxisLabelFontFile = ''
outline1Display_2.PolarAxes.LastRadialAxisTextFontFile = ''
outline1Display_2.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for mask2LUT in view renderView3
mask2LUTColorBar = GetScalarBar(mask2LUT, renderView3)
mask2LUTColorBar.Position = [0.806375442739079, 0.0747826086956522]
mask2LUTColorBar.Title = 'mask2'
mask2LUTColorBar.ComponentTitle = ''
mask2LUTColorBar.HorizontalTitle = 1
mask2LUTColorBar.TitleFontFamily = 'Courier'
mask2LUTColorBar.TitleFontFile = ''
mask2LUTColorBar.TitleBold = 1
mask2LUTColorBar.TitleFontSize = 20
mask2LUTColorBar.LabelFontFile = ''
mask2LUTColorBar.AutomaticLabelFormat = 0
mask2LUTColorBar.LabelFormat = '%-#6.1f'
mask2LUTColorBar.AddRangeLabels = 0
mask2LUTColorBar.RangeLabelFormat = '%-#6.1f'
mask2LUTColorBar.ScalarBarLength = 0.3299999999999996

# set color bar visibility
mask2LUTColorBar.Visibility = 1

# show color legend
p_vtsDisplay_2.SetScalarBarVisibility(renderView3, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView4'
# ----------------------------------------------------------------

# show data from p_vts
p_vtsDisplay_3 = Show(p_vts, renderView4)

# get color transfer function/color map for 'cl1'
cl1LUT = GetColorTransferFunction('cl1')
cl1LUT.RGBPoints = [-1.0, 0.231373, 0.298039, 0.752941, 3.5, 0.865003, 0.865003, 0.865003, 8.0, 0.705882, 0.0156863, 0.14902]
cl1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'cl1'
cl1PWF = GetOpacityTransferFunction('cl1')
cl1PWF.Points = [-1.0, 0.0, 0.5, 0.0, 8.0, 1.0, 0.5, 0.0]
cl1PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
p_vtsDisplay_3.Representation = 'Surface With Edges'
p_vtsDisplay_3.ColorArrayName = ['CELLS', 'cl1']
p_vtsDisplay_3.LookupTable = cl1LUT
p_vtsDisplay_3.Opacity = 0.5
p_vtsDisplay_3.EdgeColor = [0.3137254901960784, 0.3137254901960784, 0.3137254901960784]
p_vtsDisplay_3.OSPRayScaleFunction = 'PiecewiseFunction'
p_vtsDisplay_3.SelectOrientationVectors = 'None'
p_vtsDisplay_3.ScaleFactor = 0.1
p_vtsDisplay_3.SelectScaleArray = 'None'
p_vtsDisplay_3.GlyphType = 'Arrow'
p_vtsDisplay_3.GlyphTableIndexArray = 'None'
p_vtsDisplay_3.GaussianRadius = 0.005
p_vtsDisplay_3.SetScaleArray = [None, '']
p_vtsDisplay_3.ScaleTransferFunction = 'PiecewiseFunction'
p_vtsDisplay_3.OpacityArray = [None, '']
p_vtsDisplay_3.OpacityTransferFunction = 'PiecewiseFunction'
p_vtsDisplay_3.DataAxesGrid = 'GridAxesRepresentation'
p_vtsDisplay_3.SelectionCellLabelFontFile = ''
p_vtsDisplay_3.SelectionPointLabelFontFile = ''
p_vtsDisplay_3.PolarAxes = 'PolarAxesRepresentation'
p_vtsDisplay_3.ScalarOpacityFunction = cl1PWF
p_vtsDisplay_3.ScalarOpacityUnitDistance = 0.08839374228030177

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
p_vtsDisplay_3.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
p_vtsDisplay_3.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
p_vtsDisplay_3.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
p_vtsDisplay_3.DataAxesGrid.XTitleFontFile = ''
p_vtsDisplay_3.DataAxesGrid.YTitleFontFile = ''
p_vtsDisplay_3.DataAxesGrid.ZTitleFontFile = ''
p_vtsDisplay_3.DataAxesGrid.XLabelFontFile = ''
p_vtsDisplay_3.DataAxesGrid.YLabelFontFile = ''
p_vtsDisplay_3.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
p_vtsDisplay_3.PolarAxes.PolarAxisTitleFontFile = ''
p_vtsDisplay_3.PolarAxes.PolarAxisLabelFontFile = ''
p_vtsDisplay_3.PolarAxes.LastRadialAxisTextFontFile = ''
p_vtsDisplay_3.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from s_vtk
s_vtkDisplay_3 = Show(s_vtk, renderView4)

# trace defaults for the display properties.
s_vtkDisplay_3.Representation = 'Wireframe'
s_vtkDisplay_3.ColorArrayName = ['CELLS', 'l']
s_vtkDisplay_3.LookupTable = lLUT
s_vtkDisplay_3.PointSize = 30.0
s_vtkDisplay_3.LineWidth = 4.0
s_vtkDisplay_3.RenderPointsAsSpheres = 1
s_vtkDisplay_3.OSPRayScaleFunction = 'PiecewiseFunction'
s_vtkDisplay_3.SelectOrientationVectors = 'None'
s_vtkDisplay_3.ScaleFactor = 0.02587677538394928
s_vtkDisplay_3.SelectScaleArray = 'c'
s_vtkDisplay_3.GlyphType = 'Arrow'
s_vtkDisplay_3.GlyphTableIndexArray = 'c'
s_vtkDisplay_3.GaussianRadius = 0.001293838769197464
s_vtkDisplay_3.SetScaleArray = [None, '']
s_vtkDisplay_3.ScaleTransferFunction = 'PiecewiseFunction'
s_vtkDisplay_3.OpacityArray = [None, '']
s_vtkDisplay_3.OpacityTransferFunction = 'PiecewiseFunction'
s_vtkDisplay_3.DataAxesGrid = 'GridAxesRepresentation'
s_vtkDisplay_3.SelectionCellLabelFontFile = ''
s_vtkDisplay_3.SelectionPointLabelFontFile = ''
s_vtkDisplay_3.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_vtkDisplay_3.DataAxesGrid.XTitleFontFile = ''
s_vtkDisplay_3.DataAxesGrid.YTitleFontFile = ''
s_vtkDisplay_3.DataAxesGrid.ZTitleFontFile = ''
s_vtkDisplay_3.DataAxesGrid.XLabelFontFile = ''
s_vtkDisplay_3.DataAxesGrid.YLabelFontFile = ''
s_vtkDisplay_3.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_vtkDisplay_3.PolarAxes.PolarAxisTitleFontFile = ''
s_vtkDisplay_3.PolarAxes.PolarAxisLabelFontFile = ''
s_vtkDisplay_3.PolarAxes.LastRadialAxisTextFontFile = ''
s_vtkDisplay_3.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from outline1
outline1Display_3 = Show(outline1, renderView4)

# trace defaults for the display properties.
outline1Display_3.Representation = 'Wireframe'
outline1Display_3.ColorArrayName = [None, '']
outline1Display_3.OSPRayScaleFunction = 'PiecewiseFunction'
outline1Display_3.SelectOrientationVectors = 'None'
outline1Display_3.ScaleFactor = 0.1
outline1Display_3.SelectScaleArray = 'None'
outline1Display_3.GlyphType = 'Arrow'
outline1Display_3.GlyphTableIndexArray = 'None'
outline1Display_3.GaussianRadius = 0.005
outline1Display_3.SetScaleArray = [None, '']
outline1Display_3.ScaleTransferFunction = 'PiecewiseFunction'
outline1Display_3.OpacityArray = [None, '']
outline1Display_3.OpacityTransferFunction = 'PiecewiseFunction'
outline1Display_3.DataAxesGrid = 'GridAxesRepresentation'
outline1Display_3.SelectionCellLabelFontFile = ''
outline1Display_3.SelectionPointLabelFontFile = ''
outline1Display_3.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
outline1Display_3.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
outline1Display_3.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
outline1Display_3.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
outline1Display_3.DataAxesGrid.XTitleFontFile = ''
outline1Display_3.DataAxesGrid.YTitleFontFile = ''
outline1Display_3.DataAxesGrid.ZTitleFontFile = ''
outline1Display_3.DataAxesGrid.XLabelFontFile = ''
outline1Display_3.DataAxesGrid.YLabelFontFile = ''
outline1Display_3.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
outline1Display_3.PolarAxes.PolarAxisTitleFontFile = ''
outline1Display_3.PolarAxes.PolarAxisLabelFontFile = ''
outline1Display_3.PolarAxes.LastRadialAxisTextFontFile = ''
outline1Display_3.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for cl1LUT in view renderView4
cl1LUTColorBar = GetScalarBar(cl1LUT, renderView4)
cl1LUTColorBar.Position = [0.8087367178276268, 0.07130434782608694]
cl1LUTColorBar.Title = 'cl1'
cl1LUTColorBar.ComponentTitle = ''
cl1LUTColorBar.HorizontalTitle = 1
cl1LUTColorBar.TitleFontFamily = 'Courier'
cl1LUTColorBar.TitleFontFile = ''
cl1LUTColorBar.TitleBold = 1
cl1LUTColorBar.TitleFontSize = 20
cl1LUTColorBar.LabelFontFile = ''
cl1LUTColorBar.AutomaticLabelFormat = 0
cl1LUTColorBar.LabelFormat = '%-#6.1f'
cl1LUTColorBar.AddRangeLabels = 0
cl1LUTColorBar.ScalarBarLength = 0.33000000000000007

# set color bar visibility
cl1LUTColorBar.Visibility = 1

# show color legend
p_vtsDisplay_3.SetScalarBarVisibility(renderView4, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView5'
# ----------------------------------------------------------------

# show data from s_vtk
s_vtkDisplay_4 = Show(s_vtk, renderView5)

# trace defaults for the display properties.
s_vtkDisplay_4.Representation = 'Wireframe'
s_vtkDisplay_4.ColorArrayName = ['CELLS', 'l']
s_vtkDisplay_4.LookupTable = lLUT
s_vtkDisplay_4.Opacity = 0.99
s_vtkDisplay_4.PointSize = 30.0
s_vtkDisplay_4.LineWidth = 4.0
s_vtkDisplay_4.RenderPointsAsSpheres = 1
s_vtkDisplay_4.OSPRayScaleFunction = 'PiecewiseFunction'
s_vtkDisplay_4.SelectOrientationVectors = 'None'
s_vtkDisplay_4.ScaleFactor = 0.02587677538394928
s_vtkDisplay_4.SelectScaleArray = 'c'
s_vtkDisplay_4.GlyphType = 'Arrow'
s_vtkDisplay_4.GlyphTableIndexArray = 'c'
s_vtkDisplay_4.GaussianRadius = 0.001293838769197464
s_vtkDisplay_4.SetScaleArray = [None, '']
s_vtkDisplay_4.ScaleTransferFunction = 'PiecewiseFunction'
s_vtkDisplay_4.OpacityArray = [None, '']
s_vtkDisplay_4.OpacityTransferFunction = 'PiecewiseFunction'
s_vtkDisplay_4.DataAxesGrid = 'GridAxesRepresentation'
s_vtkDisplay_4.SelectionCellLabelFontFile = ''
s_vtkDisplay_4.SelectionPointLabelFontFile = ''
s_vtkDisplay_4.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_vtkDisplay_4.DataAxesGrid.XTitleFontFile = ''
s_vtkDisplay_4.DataAxesGrid.YTitleFontFile = ''
s_vtkDisplay_4.DataAxesGrid.ZTitleFontFile = ''
s_vtkDisplay_4.DataAxesGrid.XLabelFontFile = ''
s_vtkDisplay_4.DataAxesGrid.YLabelFontFile = ''
s_vtkDisplay_4.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_vtkDisplay_4.PolarAxes.PolarAxisTitleFontFile = ''
s_vtkDisplay_4.PolarAxes.PolarAxisLabelFontFile = ''
s_vtkDisplay_4.PolarAxes.LastRadialAxisTextFontFile = ''
s_vtkDisplay_4.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from threshold1
threshold1Display = Show(threshold1, renderView5)

# get color transfer function/color map for 'cl0'
cl0LUT = GetColorTransferFunction('cl0')
cl0LUT.RGBPoints = [-1.0, 0.0, 0.0, 1.0, 11.0, 1.0, 0.0, 0.0]
cl0LUT.ColorSpace = 'HSV'
cl0LUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
cl0LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'cl0'
cl0PWF = GetOpacityTransferFunction('cl0')
cl0PWF.Points = [-1.0, 0.0, 0.5, 0.0, 11.0, 1.0, 0.5, 0.0]
cl0PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['CELLS', 'cl0']
threshold1Display.LookupTable = cl0LUT
threshold1Display.Opacity = 0.25
threshold1Display.EdgeColor = [0.5019607843137255, 0.5019607843137255, 0.5019607843137255]
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'None'
threshold1Display.ScaleFactor = 0.1
threshold1Display.SelectScaleArray = 'None'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'None'
threshold1Display.GaussianRadius = 0.005
threshold1Display.SetScaleArray = [None, '']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = [None, '']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.SelectionCellLabelFontFile = ''
threshold1Display.SelectionPointLabelFontFile = ''
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityFunction = cl0PWF
threshold1Display.ScalarOpacityUnitDistance = 0.14034200668144686

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
threshold1Display.DataAxesGrid.XTitleFontFile = ''
threshold1Display.DataAxesGrid.YTitleFontFile = ''
threshold1Display.DataAxesGrid.ZTitleFontFile = ''
threshold1Display.DataAxesGrid.XLabelFontFile = ''
threshold1Display.DataAxesGrid.YLabelFontFile = ''
threshold1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
threshold1Display.PolarAxes.PolarAxisTitleFontFile = ''
threshold1Display.PolarAxes.PolarAxisLabelFontFile = ''
threshold1Display.PolarAxes.LastRadialAxisTextFontFile = ''
threshold1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from extractEdges1
extractEdges1Display = Show(extractEdges1, renderView5)

# trace defaults for the display properties.
extractEdges1Display.Representation = 'Wireframe'
extractEdges1Display.AmbientColor = [0.5019607843137255, 0.5019607843137255, 0.5019607843137255]
extractEdges1Display.ColorArrayName = ['POINTS', '']
extractEdges1Display.Opacity = 0.2
extractEdges1Display.PointSize = 30.0
extractEdges1Display.LineWidth = 3.0
extractEdges1Display.RenderPointsAsSpheres = 1
extractEdges1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractEdges1Display.SelectOrientationVectors = 'None'
extractEdges1Display.ScaleFactor = 0.1
extractEdges1Display.SelectScaleArray = 'None'
extractEdges1Display.GlyphType = 'Arrow'
extractEdges1Display.GlyphTableIndexArray = 'None'
extractEdges1Display.GaussianRadius = 0.005
extractEdges1Display.SetScaleArray = [None, '']
extractEdges1Display.ScaleTransferFunction = 'PiecewiseFunction'
extractEdges1Display.OpacityArray = [None, '']
extractEdges1Display.OpacityTransferFunction = 'PiecewiseFunction'
extractEdges1Display.DataAxesGrid = 'GridAxesRepresentation'
extractEdges1Display.SelectionCellLabelFontFile = ''
extractEdges1Display.SelectionPointLabelFontFile = ''
extractEdges1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
extractEdges1Display.DataAxesGrid.XTitleFontFile = ''
extractEdges1Display.DataAxesGrid.YTitleFontFile = ''
extractEdges1Display.DataAxesGrid.ZTitleFontFile = ''
extractEdges1Display.DataAxesGrid.XLabelFontFile = ''
extractEdges1Display.DataAxesGrid.YLabelFontFile = ''
extractEdges1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
extractEdges1Display.PolarAxes.PolarAxisTitleFontFile = ''
extractEdges1Display.PolarAxes.PolarAxisLabelFontFile = ''
extractEdges1Display.PolarAxes.LastRadialAxisTextFontFile = ''
extractEdges1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'l'
lPWF = GetOpacityTransferFunction('l')
lPWF.Points = [0.0, 1.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]
lPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(threshold1)
# ----------------------------------------------------------------