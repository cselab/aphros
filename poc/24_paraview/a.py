# state file generated using paraview version 5.7.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.7.0
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
renderView1.ViewSize = [2158, 1180]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5, 0.39997434616088867, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.5, 0.39997434616088867, 3.638838601971638]
renderView1.CameraFocalPoint = [0.5, 0.39997434616088867, 0.5]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.8123912096932305
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

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
sm_ = LegacyVTKReader(FileNames=['/home/kpetr/s/euler_scratch/ch/pasc/rise/nx064blayer8cflst01/sm_0000.vtk', '/home/kpetr/s/euler_scratch/ch/pasc/rise/nx064blayer8cflst01/sm_0001.vtk'])

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from sm_
sm_Display = Show(sm_, renderView1)

# get color transfer function/color map for 'c'
cLUT = GetColorTransferFunction('c')
cLUT.RGBPoints = [2004.0, 0.231373, 0.298039, 0.752941, 31526534.000000007, 0.865003, 0.865003, 0.865003, 63051064.0, 0.705882, 0.0156863, 0.14902]
cLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
sm_Display.Representation = 'Surface'
sm_Display.AmbientColor = [0.0, 0.0, 0.0]
sm_Display.ColorArrayName = ['CELLS', 'c']
sm_Display.LookupTable = cLUT
sm_Display.OSPRayScaleArray = 'nn'
sm_Display.OSPRayScaleFunction = 'PiecewiseFunction'
sm_Display.SelectOrientationVectors = 'nn'
sm_Display.ScaleFactor = 0.1
sm_Display.SelectScaleArray = 'c'
sm_Display.GlyphType = 'Arrow'
sm_Display.GlyphTableIndexArray = 'c'
sm_Display.GaussianRadius = 0.005
sm_Display.SetScaleArray = ['POINTS', 'nn']
sm_Display.ScaleTransferFunction = 'PiecewiseFunction'
sm_Display.OpacityArray = ['POINTS', 'nn']
sm_Display.OpacityTransferFunction = 'PiecewiseFunction'
sm_Display.DataAxesGrid = 'GridAxesRepresentation'
sm_Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
sm_Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sm_Display.ScaleTransferFunction.Points = [-0.9512056112289429, 1.0, 0.5, 0.0, 0.9429839849472046, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sm_Display.OpacityTransferFunction.Points = [-0.9512056112289429, 1.0, 0.5, 0.0, 0.9429839849472046, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
sm_Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
sm_Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
sm_Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
sm_Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
sm_Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
sm_Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
sm_Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
sm_Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
sm_Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
sm_Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
sm_Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for cLUT in view renderView1
cLUTColorBar = GetScalarBar(cLUT, renderView1)
cLUTColorBar.Title = 'c'
cLUTColorBar.ComponentTitle = ''
cLUTColorBar.HorizontalTitle = 1
cLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
cLUTColorBar.TitleFontFamily = 'Courier'
cLUTColorBar.TitleBold = 1
cLUTColorBar.TitleFontSize = 20
cLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
cLUTColorBar.AutomaticLabelFormat = 0
cLUTColorBar.LabelFormat = '%-#6.1f'
cLUTColorBar.AddRangeLabels = 0
cLUTColorBar.ScalarBarLength = 0.33000000000000007

# set color bar visibility
cLUTColorBar.Visibility = 1

# show color legend
sm_Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'c'
cPWF = GetOpacityTransferFunction('c')
cPWF.Points = [2004.0, 1.0, 0.5, 0.0, 63051064.0, 1.0, 0.5, 0.0]
cPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(sm_)
# ----------------------------------------------------------------
