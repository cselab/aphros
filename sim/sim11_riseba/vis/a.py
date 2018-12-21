# state file generated using paraview version 5.6.0-RC2

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.6.0-RC2
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
renderView1.ViewSize = [640, 360]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [1.0592000070190428, 0.16607749462127686, 0.0078125]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.9886528578892481, 0.23633099725623352, 1.1266127606454457]
renderView1.CameraFocalPoint = [0.9886528578892481, 0.23633099725623352, 0.0078125]
renderView1.CameraParallelScale = 0.23931141745512827
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
c2g3l4stxt = CSVReader(FileName=['/home/petr/s/euler_sim10_riseba/c2g3l4s.txt'])
c2g3l4stxt.HaveHeaders = 0
c2g3l4stxt.FieldDelimiterCharacters = ' '
c2g3l4stxt.AddTabFieldDelimiter = 1
c2g3l4stxt.MergeConsecutiveDelimiters = 1

# create a new 'CSV Reader'
rising2ref = CSVReader(FileName=['/home/petr/s/euler_sim10_riseba/rising2.ref'])
rising2ref.HaveHeaders = 0
rising2ref.FieldDelimiterCharacters = ' '

# create a new 'CSV Reader'
c1g3l4stxt = CSVReader(FileName=['/home/petr/s/euler_sim10_riseba/c1g3l4s.txt'])
c1g3l4stxt.HaveHeaders = 0
c1g3l4stxt.FieldDelimiterCharacters = ' '
c1g3l4stxt.AddTabFieldDelimiter = 1
c1g3l4stxt.MergeConsecutiveDelimiters = 1

# create a new 'Legacy VTK Reader'
c1ch = LegacyVTKReader(FileNames=['/home/petr/s/euler_sim10_riseba/s_0000.vtk', '/home/petr/s/euler_sim10_riseba/s_0001.vtk', '/home/petr/s/euler_sim10_riseba/s_0002.vtk', '/home/petr/s/euler_sim10_riseba/s_0003.vtk', '/home/petr/s/euler_sim10_riseba/s_0004.vtk', '/home/petr/s/euler_sim10_riseba/s_0005.vtk', '/home/petr/s/euler_sim10_riseba/s_0006.vtk', '/home/petr/s/euler_sim10_riseba/s_0007.vtk', '/home/petr/s/euler_sim10_riseba/s_0008.vtk', '/home/petr/s/euler_sim10_riseba/s_0009.vtk', '/home/petr/s/euler_sim10_riseba/s_0010.vtk', '/home/petr/s/euler_sim10_riseba/s_0011.vtk', '/home/petr/s/euler_sim10_riseba/s_0012.vtk', '/home/petr/s/euler_sim10_riseba/s_0013.vtk', '/home/petr/s/euler_sim10_riseba/s_0014.vtk', '/home/petr/s/euler_sim10_riseba/s_0015.vtk', '/home/petr/s/euler_sim10_riseba/s_0016.vtk', '/home/petr/s/euler_sim10_riseba/s_0017.vtk', '/home/petr/s/euler_sim10_riseba/s_0018.vtk', '/home/petr/s/euler_sim10_riseba/s_0019.vtk', '/home/petr/s/euler_sim10_riseba/s_0020.vtk', '/home/petr/s/euler_sim10_riseba/s_0021.vtk', '/home/petr/s/euler_sim10_riseba/s_0022.vtk', '/home/petr/s/euler_sim10_riseba/s_0023.vtk', '/home/petr/s/euler_sim10_riseba/s_0024.vtk', '/home/petr/s/euler_sim10_riseba/s_0025.vtk', '/home/petr/s/euler_sim10_riseba/s_0026.vtk', '/home/petr/s/euler_sim10_riseba/s_0027.vtk', '/home/petr/s/euler_sim10_riseba/s_0028.vtk', '/home/petr/s/euler_sim10_riseba/s_0029.vtk', '/home/petr/s/euler_sim10_riseba/s_0030.vtk'])

# create a new 'Table To Points'
tableToPoints4 = TableToPoints(Input=c2g3l4stxt)
tableToPoints4.XColumn = 'Field 1'
tableToPoints4.YColumn = 'Field 0'
tableToPoints4.ZColumn = 'Field 0'

# create a new 'CSV Reader'
risingref = CSVReader(FileName=['/home/petr/s/ch/sim/sim10_riseba/rising.ref'])
risingref.HaveHeaders = 0
risingref.FieldDelimiterCharacters = ' '

# create a new 'Table To Points'
c1ba = TableToPoints(Input=risingref)
c1ba.XColumn = 'Field 0'
c1ba.YColumn = 'Field 1'
c1ba.ZColumn = 'Field 0'
c1ba.a2DPoints = 1

# create a new 'Legacy VTK Reader'
c2ch = LegacyVTKReader(FileNames=['/home/petr/s/euler_sim10_riseba/c1/s_0000.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0001.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0002.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0003.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0004.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0005.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0006.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0007.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0008.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0009.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0010.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0011.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0012.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0013.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0014.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0015.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0016.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0017.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0018.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0019.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0020.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0021.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0022.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0023.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0024.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0025.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0026.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0027.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0028.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0029.vtk', '/home/petr/s/euler_sim10_riseba/c1/s_0030.vtk'])

# create a new 'Table To Points'
c2ba = TableToPoints(Input=rising2ref)
c2ba.XColumn = 'Field 0'
c2ba.YColumn = 'Field 1'
c2ba.ZColumn = 'Field 0'
c2ba.a2DPoints = 1

# create a new 'Table To Points'
tableToPoints3 = TableToPoints(Input=c1g3l4stxt)
tableToPoints3.XColumn = 'Field 1'
tableToPoints3.YColumn = 'Field 0'
tableToPoints3.ZColumn = 'Field 0'
tableToPoints3.a2DPoints = 1

# create a new 'Transform'
c1mn = Transform(Input=tableToPoints3)
c1mn.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
c1mn.Transform.Translate = [0.0, -0.5, 0.0]

# create a new 'Transform'
c2mn = Transform(Input=tableToPoints4)
c2mn.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
c2mn.Transform.Translate = [0.0, -0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from c1ba
c1baDisplay = Show(c1ba, renderView1)

# trace defaults for the display properties.
c1baDisplay.Representation = 'Points'
c1baDisplay.AmbientColor = [0.0, 1.0, 0.0]
c1baDisplay.ColorArrayName = [None, '']
c1baDisplay.PointSize = 5.0
c1baDisplay.LineWidth = 3.0
c1baDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
c1baDisplay.SelectOrientationVectors = 'None'
c1baDisplay.ScaleFactor = 0.07650390000000001
c1baDisplay.SelectScaleArray = 'None'
c1baDisplay.GlyphType = 'Arrow'
c1baDisplay.GlyphTableIndexArray = 'None'
c1baDisplay.GaussianRadius = 0.0038251950000000004
c1baDisplay.SetScaleArray = [None, '']
c1baDisplay.ScaleTransferFunction = 'PiecewiseFunction'
c1baDisplay.OpacityArray = [None, '']
c1baDisplay.OpacityTransferFunction = 'PiecewiseFunction'
c1baDisplay.DataAxesGrid = 'GridAxesRepresentation'
c1baDisplay.SelectionCellLabelFontFile = ''
c1baDisplay.SelectionPointLabelFontFile = ''
c1baDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
c1baDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
c1baDisplay.DataAxesGrid.XTitleFontFile = ''
c1baDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
c1baDisplay.DataAxesGrid.YTitleFontFile = ''
c1baDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
c1baDisplay.DataAxesGrid.ZTitleFontFile = ''
c1baDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
c1baDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
c1baDisplay.DataAxesGrid.XLabelFontFile = ''
c1baDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
c1baDisplay.DataAxesGrid.YLabelFontFile = ''
c1baDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
c1baDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
c1baDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
c1baDisplay.PolarAxes.PolarAxisTitleFontFile = ''
c1baDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
c1baDisplay.PolarAxes.PolarAxisLabelFontFile = ''
c1baDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
c1baDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
c1baDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
c1baDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from c1ch
c1chDisplay = Show(c1ch, renderView1)

# trace defaults for the display properties.
c1chDisplay.Representation = 'Wireframe'
c1chDisplay.AmbientColor = [0.3333333333333333, 0.3333333333333333, 1.0]
c1chDisplay.ColorArrayName = ['POINTS', '']
c1chDisplay.PointSize = 30.0
c1chDisplay.LineWidth = 3.0
c1chDisplay.RenderPointsAsSpheres = 1
c1chDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
c1chDisplay.SelectOrientationVectors = 'None'
c1chDisplay.ScaleFactor = 0.03692150115966797
c1chDisplay.SelectScaleArray = 'c'
c1chDisplay.GlyphType = 'Arrow'
c1chDisplay.GlyphTableIndexArray = 'c'
c1chDisplay.GaussianRadius = 0.0018460750579833984
c1chDisplay.SetScaleArray = [None, '']
c1chDisplay.ScaleTransferFunction = 'PiecewiseFunction'
c1chDisplay.OpacityArray = [None, '']
c1chDisplay.OpacityTransferFunction = 'PiecewiseFunction'
c1chDisplay.DataAxesGrid = 'GridAxesRepresentation'
c1chDisplay.SelectionCellLabelFontFile = ''
c1chDisplay.SelectionPointLabelFontFile = ''
c1chDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
c1chDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
c1chDisplay.DataAxesGrid.XTitleFontFile = ''
c1chDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
c1chDisplay.DataAxesGrid.YTitleFontFile = ''
c1chDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
c1chDisplay.DataAxesGrid.ZTitleFontFile = ''
c1chDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
c1chDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
c1chDisplay.DataAxesGrid.XLabelFontFile = ''
c1chDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
c1chDisplay.DataAxesGrid.YLabelFontFile = ''
c1chDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
c1chDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
c1chDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
c1chDisplay.PolarAxes.PolarAxisTitleFontFile = ''
c1chDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
c1chDisplay.PolarAxes.PolarAxisLabelFontFile = ''
c1chDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
c1chDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
c1chDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
c1chDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from c1mn
c1mnDisplay = Show(c1mn, renderView1)

# trace defaults for the display properties.
c1mnDisplay.Representation = 'Surface'
c1mnDisplay.AmbientColor = [0.3333333333333333, 0.3333333333333333, 1.0]
c1mnDisplay.ColorArrayName = [None, '']
c1mnDisplay.DiffuseColor = [1.0, 0.0, 0.0]
c1mnDisplay.PointSize = 3.0
c1mnDisplay.LineWidth = 3.0
c1mnDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
c1mnDisplay.SelectOrientationVectors = 'None'
c1mnDisplay.ScaleFactor = 0.0706143
c1mnDisplay.SelectScaleArray = 'None'
c1mnDisplay.GlyphType = 'Arrow'
c1mnDisplay.GlyphTableIndexArray = 'None'
c1mnDisplay.GaussianRadius = 0.0035307150000000002
c1mnDisplay.SetScaleArray = [None, '']
c1mnDisplay.ScaleTransferFunction = 'PiecewiseFunction'
c1mnDisplay.OpacityArray = [None, '']
c1mnDisplay.OpacityTransferFunction = 'PiecewiseFunction'
c1mnDisplay.DataAxesGrid = 'GridAxesRepresentation'
c1mnDisplay.SelectionCellLabelFontFile = ''
c1mnDisplay.SelectionPointLabelFontFile = ''
c1mnDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
c1mnDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
c1mnDisplay.DataAxesGrid.XTitleFontFile = ''
c1mnDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
c1mnDisplay.DataAxesGrid.YTitleFontFile = ''
c1mnDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
c1mnDisplay.DataAxesGrid.ZTitleFontFile = ''
c1mnDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
c1mnDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
c1mnDisplay.DataAxesGrid.XLabelFontFile = ''
c1mnDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
c1mnDisplay.DataAxesGrid.YLabelFontFile = ''
c1mnDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
c1mnDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
c1mnDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
c1mnDisplay.PolarAxes.PolarAxisTitleFontFile = ''
c1mnDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
c1mnDisplay.PolarAxes.PolarAxisLabelFontFile = ''
c1mnDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
c1mnDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
c1mnDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
c1mnDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(c1mn)
# ----------------------------------------------------------------