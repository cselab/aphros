# state file generated using paraview version 5.5.2

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1376, 972]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [-0.077902768708841, 0.9194849591490677, 0.9356836987833435]
renderView1.CameraFocalPoint = [0.5000000000000003, 0.4999999999999995, 0.49999999999999994]
renderView1.CameraViewUp = [-0.6669417150515764, -0.16387759336564117, -0.726865106547802]
renderView1.CameraParallelScale = 0.21650635094610965
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
partit_0000csv = CSVReader(FileName=['/home/kpetr/s/ch/sim/partstr/iter/partit_0000.csv'])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=partit_0000csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Legacy VTK Reader'
s_0000vtk = LegacyVTKReader(FileNames=['/home/kpetr/s/ch/sim/partstr/iter/s_0000.vtk'])

# create a new 'Glyph'
glyph1 = Glyph(Input=tableToPoints1,
    GlyphType='Sphere')
glyph1.Scalars = ['POINTS', 'None']
glyph1.Vectors = ['POINTS', 'None']
glyph1.ScaleFactor = 0.01
glyph1.GlyphTransform = 'Transform2'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from s_0000vtk
s_0000vtkDisplay = Show(s_0000vtk, renderView1)

# trace defaults for the display properties.
s_0000vtkDisplay.Representation = 'Surface'
s_0000vtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.ColorArrayName = ['POINTS', '']
s_0000vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s_0000vtkDisplay.SelectOrientationVectors = 'None'
s_0000vtkDisplay.ScaleFactor = 0.025
s_0000vtkDisplay.SelectScaleArray = 'c'
s_0000vtkDisplay.GlyphType = 'Arrow'
s_0000vtkDisplay.GlyphTableIndexArray = 'c'
s_0000vtkDisplay.GaussianRadius = 0.00125
s_0000vtkDisplay.SetScaleArray = [None, '']
s_0000vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
s_0000vtkDisplay.OpacityArray = [None, '']
s_0000vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
s_0000vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
s_0000vtkDisplay.SelectionCellLabelFontFile = ''
s_0000vtkDisplay.SelectionPointLabelFontFile = ''
s_0000vtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
s_0000vtkDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
s_0000vtkDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
s_0000vtkDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_0000vtkDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.DataAxesGrid.XTitleFontFile = ''
s_0000vtkDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.DataAxesGrid.YTitleFontFile = ''
s_0000vtkDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.DataAxesGrid.ZTitleFontFile = ''
s_0000vtkDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.DataAxesGrid.XLabelFontFile = ''
s_0000vtkDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.DataAxesGrid.YLabelFontFile = ''
s_0000vtkDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_0000vtkDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
s_0000vtkDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
s_0000vtkDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
s_0000vtkDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
s_0000vtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from tableToPoints1
tableToPoints1Display = Show(tableToPoints1, renderView1)

# trace defaults for the display properties.
tableToPoints1Display.Representation = 'Surface'
tableToPoints1Display.ColorArrayName = [None, '']
tableToPoints1Display.OSPRayScaleArray = 'c'
tableToPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display.SelectOrientationVectors = 'None'
tableToPoints1Display.ScaleFactor = 0.024892455074563547
tableToPoints1Display.SelectScaleArray = 'None'
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.GlyphTableIndexArray = 'None'
tableToPoints1Display.GaussianRadius = 0.0012446227537281774
tableToPoints1Display.SetScaleArray = ['POINTS', 'c']
tableToPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.OpacityArray = ['POINTS', 'c']
tableToPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display.SelectionCellLabelFontFile = ''
tableToPoints1Display.SelectionPointLabelFontFile = ''
tableToPoints1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
tableToPoints1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tableToPoints1Display.ScaleTransferFunction.Points = [6006006.0, 1.0, 0.5, 0.0, 9009009.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tableToPoints1Display.OpacityTransferFunction.Points = [6006006.0, 1.0, 0.5, 0.0, 9009009.0, 1.0, 0.5, 0.0]

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

# show data from glyph1
glyph1Display = Show(glyph1, renderView1)

# get color transfer function/color map for 'c'
cLUT = GetColorTransferFunction('c')
cLUT.RGBPoints = [6006006.0, 0.231373, 0.298039, 0.752941, 7507508.0, 0.865003, 0.865003, 0.865003, 9009010.0, 0.705882, 0.0156863, 0.14902]
cLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'c']
glyph1Display.LookupTable = cLUT
glyph1Display.OSPRayScaleArray = 'Normals'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 0.02738170325756073
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.0013690851628780365
glyph1Display.SetScaleArray = ['POINTS', 'Normals']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Normals']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.SelectionCellLabelFontFile = ''
glyph1Display.SelectionPointLabelFontFile = ''
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-0.9749279618263245, 1.0, 0.5, 0.0, 0.9749279618263245, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-0.9749279618263245, 1.0, 0.5, 0.0, 0.9749279618263245, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleFontFile = ''
glyph1Display.DataAxesGrid.YTitleFontFile = ''
glyph1Display.DataAxesGrid.ZTitleFontFile = ''
glyph1Display.DataAxesGrid.XLabelFontFile = ''
glyph1Display.DataAxesGrid.YLabelFontFile = ''
glyph1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleFontFile = ''
glyph1Display.PolarAxes.PolarAxisLabelFontFile = ''
glyph1Display.PolarAxes.LastRadialAxisTextFontFile = ''
glyph1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for cLUT in view renderView1
cLUTColorBar = GetScalarBar(cLUT, renderView1)
cLUTColorBar.Title = 'c'
cLUTColorBar.ComponentTitle = ''
cLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
cLUTColorBar.TitleFontFile = ''
cLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
cLUTColorBar.LabelFontFile = ''

# set color bar visibility
cLUTColorBar.Visibility = 1

# show color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'c'
cPWF = GetOpacityTransferFunction('c')
cPWF.Points = [6006006.0, 1.0, 0.5, 0.0, 9009010.0, 1.0, 0.5, 0.0]
cPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(s_0000vtk)
# ----------------------------------------------------------------