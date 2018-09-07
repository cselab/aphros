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
renderView1.ViewSize = [2171, 918]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5, 0.5, 0.0]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.5, 0.5, 2.7320508075688776]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476
renderView1.CameraParallelProjection = 1
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

# create a new 'Legacy VTK Reader'
s_0 = LegacyVTKReader(FileNames=['/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0000.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0001.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0002.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0003.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0004.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0005.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0006.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0007.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0008.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0009.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0010.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0011.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0012.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0013.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0014.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0015.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0016.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0017.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0018.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0019.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0020.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0021.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0022.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0023.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0024.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0025.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0026.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0027.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0028.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0029.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0030.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0031.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0032.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0033.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0034.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0035.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0036.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0037.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0038.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0039.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0040.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0041.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0042.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0043.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0044.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0045.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0046.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0047.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0048.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0049.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0050.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0051.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0052.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0053.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0054.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0055.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0056.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0057.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0058.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0059.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0060.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0061.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0062.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0063.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0064.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0065.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0066.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0067.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0068.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0069.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0070.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0071.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0072.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0073.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0074.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0075.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0076.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0077.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0078.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0079.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0080.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0081.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0082.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0083.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0084.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0085.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0086.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0087.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0088.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0089.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0090.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0091.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0092.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0093.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0094.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0095.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0096.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0097.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0098.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0099.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ch/s_0100.vtk'])

# create a new 'Legacy VTK Reader'
u_0 = LegacyVTKReader(FileNames=['/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0000.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0001.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0002.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0003.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0004.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0005.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0006.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0007.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0008.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0009.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0010.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0011.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0012.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0013.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0014.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0015.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0016.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0017.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0018.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0019.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0020.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0021.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0022.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0023.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0024.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0025.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0026.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0027.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0028.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0029.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0030.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0031.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0032.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0033.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0034.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0035.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0036.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0037.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0038.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0039.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0040.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0041.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0042.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0043.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0044.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0045.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0046.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0047.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0048.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0049.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0050.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0051.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0052.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0053.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0054.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0055.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0056.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0057.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0058.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0059.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0060.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0061.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0062.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0063.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0064.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0065.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0066.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0067.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0068.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0069.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0070.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0071.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0072.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0073.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0074.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0075.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0076.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0077.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0078.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0079.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0080.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0081.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0082.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0083.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0084.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0085.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0086.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0087.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0088.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0089.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0090.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0091.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0092.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0093.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0094.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0095.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0096.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0097.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0098.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0099.vtk', '/home/kpetr/s/ch/sim/sim10_univel3d/ge/u_0100.vtk'])

# create a new 'Transform'
transform1 = Transform(Input=u_0)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [0.5, 0.5, 0.0]

# create a new 'Contour'
contour1 = Contour(Input=transform1)
contour1.ContourBy = ['POINTS', 'T']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Outline'
outline1 = Outline(Input=transform1)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from s_0
s_0Display = Show(s_0, renderView1)

# trace defaults for the display properties.
s_0Display.Representation = 'Wireframe'
s_0Display.AmbientColor = [0.12156862745098039, 0.4666666666666667, 0.7058823529411765]
s_0Display.ColorArrayName = ['POINTS', '']
s_0Display.LineWidth = 3.0
s_0Display.OSPRayScaleFunction = 'PiecewiseFunction'
s_0Display.SelectOrientationVectors = 'None'
s_0Display.ScaleFactor = 0.06021709740161896
s_0Display.SelectScaleArray = 'c'
s_0Display.GlyphType = 'Arrow'
s_0Display.GlyphTableIndexArray = 'c'
s_0Display.GaussianRadius = 0.003010854870080948
s_0Display.SetScaleArray = [None, '']
s_0Display.ScaleTransferFunction = 'PiecewiseFunction'
s_0Display.OpacityArray = [None, '']
s_0Display.OpacityTransferFunction = 'PiecewiseFunction'
s_0Display.DataAxesGrid = 'GridAxesRepresentation'
s_0Display.SelectionCellLabelFontFile = ''
s_0Display.SelectionPointLabelFontFile = ''
s_0Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_0Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
s_0Display.DataAxesGrid.XTitleFontFile = ''
s_0Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
s_0Display.DataAxesGrid.YTitleFontFile = ''
s_0Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
s_0Display.DataAxesGrid.ZTitleFontFile = ''
s_0Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
s_0Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
s_0Display.DataAxesGrid.XLabelFontFile = ''
s_0Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
s_0Display.DataAxesGrid.YLabelFontFile = ''
s_0Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
s_0Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_0Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
s_0Display.PolarAxes.PolarAxisTitleFontFile = ''
s_0Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
s_0Display.PolarAxes.PolarAxisLabelFontFile = ''
s_0Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
s_0Display.PolarAxes.LastRadialAxisTextFontFile = ''
s_0Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
s_0Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from transform1
transform1Display = Show(transform1, renderView1)

# trace defaults for the display properties.
transform1Display.Representation = 'Wireframe'
transform1Display.AmbientColor = [0.0, 0.0, 0.0]
transform1Display.ColorArrayName = ['POINTS', '']
transform1Display.Opacity = 0.25
transform1Display.OSPRayScaleArray = 'P'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 0.1
transform1Display.SelectScaleArray = 'None'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'None'
transform1Display.GaussianRadius = 0.005
transform1Display.SetScaleArray = ['POINTS', 'P']
transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display.OpacityArray = ['POINTS', 'P']
transform1Display.OpacityTransferFunction = 'PiecewiseFunction'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.SelectionCellLabelFontFile = ''
transform1Display.SelectionPointLabelFontFile = ''
transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.ScalarOpacityUnitDistance = 0.22272467953508485

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-0.161638, 0.0, 0.5, 0.0, 0.6202350000000001, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-0.161638, 0.0, 0.5, 0.0, 0.6202350000000001, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
transform1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
transform1Display.DataAxesGrid.XTitleFontFile = ''
transform1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
transform1Display.DataAxesGrid.YTitleFontFile = ''
transform1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
transform1Display.DataAxesGrid.ZTitleFontFile = ''
transform1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
transform1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
transform1Display.DataAxesGrid.XLabelFontFile = ''
transform1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
transform1Display.DataAxesGrid.YLabelFontFile = ''
transform1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
transform1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
transform1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
transform1Display.PolarAxes.PolarAxisTitleFontFile = ''
transform1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
transform1Display.PolarAxes.PolarAxisLabelFontFile = ''
transform1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
transform1Display.PolarAxes.LastRadialAxisTextFontFile = ''
transform1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
transform1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.ColorArrayName = ['POINTS', '']
contour1Display.DiffuseColor = [1.0, 0.4980392156862745, 0.054901960784313725]
contour1Display.LineWidth = 3.0
contour1Display.OSPRayScaleArray = 'P'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.058940085964912285
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.GaussianRadius = 0.002947004298245614
contour1Display.SetScaleArray = ['POINTS', 'P']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'P']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.SelectionCellLabelFontFile = ''
contour1Display.SelectionPointLabelFontFile = ''
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [0.1884364335609578, 0.0, 0.5, 0.0, 0.2774627792019816, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.1884364335609578, 0.0, 0.5, 0.0, 0.2774627792019816, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
contour1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
contour1Display.DataAxesGrid.XTitleFontFile = ''
contour1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
contour1Display.DataAxesGrid.YTitleFontFile = ''
contour1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
contour1Display.DataAxesGrid.ZTitleFontFile = ''
contour1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
contour1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
contour1Display.DataAxesGrid.XLabelFontFile = ''
contour1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
contour1Display.DataAxesGrid.YLabelFontFile = ''
contour1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
contour1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
contour1Display.PolarAxes.PolarAxisTitleFontFile = ''
contour1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
contour1Display.PolarAxes.PolarAxisLabelFontFile = ''
contour1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
contour1Display.PolarAxes.LastRadialAxisTextFontFile = ''
contour1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
contour1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from outline1
outline1Display = Show(outline1, renderView1)

# trace defaults for the display properties.
outline1Display.Representation = 'Surface'
outline1Display.AmbientColor = [0.0, 0.0, 0.0]
outline1Display.ColorArrayName = [None, '']
outline1Display.DiffuseColor = [0.0, 0.0, 0.0]
outline1Display.LineWidth = 2.0
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

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
outline1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
outline1Display.DataAxesGrid.XTitleFontFile = ''
outline1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
outline1Display.DataAxesGrid.YTitleFontFile = ''
outline1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
outline1Display.DataAxesGrid.ZTitleFontFile = ''
outline1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
outline1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
outline1Display.DataAxesGrid.XLabelFontFile = ''
outline1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
outline1Display.DataAxesGrid.YLabelFontFile = ''
outline1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
outline1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
outline1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
outline1Display.PolarAxes.PolarAxisTitleFontFile = ''
outline1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
outline1Display.PolarAxes.PolarAxisLabelFontFile = ''
outline1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
outline1Display.PolarAxes.LastRadialAxisTextFontFile = ''
outline1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
outline1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(transform1)
# ----------------------------------------------------------------
