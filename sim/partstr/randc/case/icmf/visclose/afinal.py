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
renderView1.ViewSize = [1000, 1000]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5409649312496185, 0.4211595505475998, 0.5226505547761917]
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.StereoType = 0
renderView1.CameraPosition = [0.24839119124607448, -0.10623063113324308, 0.3704770501303224]
renderView1.CameraFocalPoint = [0.7212417088961819, 0.6928833886588502, 0.5883919766652376]
renderView1.CameraViewUp = [0.7866517756373992, -0.5447244857787051, 0.2906100797970575]
renderView1.CameraParallelScale = 0.24685119903525776
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableOSPRay = 1
renderView1.AmbientSamples = 1
renderView1.SamplesPerPixel = 5
renderView1.ProgressivePasses = 10
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
ovtk = LegacyVTKReader(FileNames=['/home/kpetr/s/electrochem/conf/icmf2019/slides/figloc/dim3_pltr010b000/o.vtk'])

# create a new 'Programmable Filter'
prog_partcon = ProgrammableFilter(Input=ovtk)
prog_partcon.Script = """execfile("/home/kpetr/s/electrochem/ch/sim/partstr/randc/case/icmf/vis/prog_part.py")
"""
prog_partcon.RequestInformationScript = ''
prog_partcon.RequestUpdateExtentScript = ''
prog_partcon.PythonPath = ''

# create a new 'Clip'
clip_partcon = Clip(Input=prog_partcon)
clip_partcon.ClipType = 'Scalar'
clip_partcon.Scalars = ['CELLS', 'sc']
clip_partcon.Value = 3713713.24
clip_partcon.Invert = 0

# create a new 'CSV Reader'
partit_csv = CSVReader(FileName=['/home/kpetr/s/electrochem/conf/icmf2019/slides/figloc/dim3_pltr010b000/partit_0000.csv'])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=partit_csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Programmable Filter'
prog_part = ProgrammableFilter(Input=tableToPoints1)
prog_part.Script = """execfile("/home/kpetr/s/electrochem/ch/sim/partstr/randc/case/icmf/vis/prog_part.py")
"""
prog_part.RequestInformationScript = ''
prog_part.RequestUpdateExtentScript = ''
prog_part.PythonPath = ''

# create a new 'Programmable Filter'
slice2 = ProgrammableFilter(Input=prog_part)
slice2.OutputDataSetType = 'vtkPolyData'
slice2.Script = """import sys
sys.argv = [0.15, 0.3, 0.06, 9, "poly"] 
execfile("/home/kpetr/s/electrochem/ch/sim/partstr/randc/case/icmf/vis/prog_polygon.py")  """
slice2.RequestInformationScript = ''
slice2.RequestUpdateExtentScript = ''
slice2.PythonPath = ''

# create a new 'Programmable Filter'
line = ProgrammableFilter(Input=prog_part)
line.Script = """import sys
sys.argv = [0.15, 0.3, 0.06, 0, "line"] 
execfile("/home/kpetr/s/electrochem/ch/sim/partstr/randc/case/icmf/vis/prog_polygon.py")  """
line.RequestInformationScript = ''
line.RequestUpdateExtentScript = ''
line.PythonPath = ''

# create a new 'Extract Edges'
extractEdges2 = ExtractEdges(Input=slice2)

# create a new 'Legacy VTK Reader'
sp_vtk = LegacyVTKReader(FileNames=['/home/kpetr/s/electrochem/conf/icmf2019/slides/figloc/dim3_pltr010b000/sp_0000.vtk'])

# create a new 'Programmable Filter'
prog_line = ProgrammableFilter(Input=sp_vtk)
prog_line.Script = """execfile("/home/kpetr/s/electrochem/ch/sim/partstr/randc/case/icmf/vis/prog_part.py")
"""
prog_line.RequestInformationScript = ''
prog_line.RequestUpdateExtentScript = ''
prog_line.PythonPath = ''

# create a new 'Clip'
clip_line = Clip(Input=prog_line)
clip_line.ClipType = 'Scalar'
clip_line.Scalars = ['CELLS', 'sc']
clip_line.Value = 7807807.8
clip_line.Invert = 0

# create a new 'Programmable Filter'
slice = ProgrammableFilter(Input=prog_part)
slice.OutputDataSetType = 'vtkPolyData'
slice.Script = """import sys
sys.argv = [0.15, 0.3, 0.06, 0, "poly"] 
execfile("/home/kpetr/s/electrochem/ch/sim/partstr/randc/case/icmf/vis/prog_polygon.py")  """
slice.RequestInformationScript = ''
slice.RequestUpdateExtentScript = ''
slice.PythonPath = ''

# create a new 'Extract Edges'
extractEdges1 = ExtractEdges(Input=slice)

# create a new 'Legacy VTK Reader'
s_vtk = LegacyVTKReader(FileNames=['/home/kpetr/s/electrochem/conf/icmf2019/slides/figloc/dim3_pltr010b000/s_0000.vtk'])

# create a new 'Clip'
clip1 = Clip(Input=s_vtk)
clip1.ClipType = 'Plane'
clip1.Scalars = ['CELLS', 'c']
clip1.Value = 8007007.5
clip1.Crinkleclip = 1

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.5788521482325462, 0.408296197740581, 0.5227404940354464]
clip1.ClipType.Normal = [-0.4573763428286224, 0.8819643456684312, -0.11377949723201745]

# create a new 'Clip'
clippart = Clip(Input=prog_part)
clippart.ClipType = 'Scalar'
clippart.Scalars = ['POINTS', 'sp']
clippart.Value = 4344343.96
clippart.Invert = 0

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from clip_partcon
clip_partconDisplay = Show(clip_partcon, renderView1)

# trace defaults for the display properties.
clip_partconDisplay.Representation = 'Surface'
clip_partconDisplay.AmbientColor = [1.0, 0.0, 0.0]
clip_partconDisplay.ColorArrayName = ['POINTS', '']
clip_partconDisplay.DiffuseColor = [1.0, 0.4980392156862745, 0.054901960784313725]
clip_partconDisplay.LineWidth = 2.0
clip_partconDisplay.RenderLinesAsTubes = 1
clip_partconDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clip_partconDisplay.SelectOrientationVectors = 'None'
clip_partconDisplay.ScaleFactor = 0.01659708619117737
clip_partconDisplay.SelectScaleArray = 'None'
clip_partconDisplay.GlyphType = 'Arrow'
clip_partconDisplay.GlyphTableIndexArray = 'None'
clip_partconDisplay.GaussianRadius = 0.0008298543095588684
clip_partconDisplay.SetScaleArray = [None, '']
clip_partconDisplay.ScaleTransferFunction = 'PiecewiseFunction'
clip_partconDisplay.OpacityArray = [None, '']
clip_partconDisplay.OpacityTransferFunction = 'PiecewiseFunction'
clip_partconDisplay.DataAxesGrid = 'GridAxesRepresentation'
clip_partconDisplay.SelectionCellLabelFontFile = ''
clip_partconDisplay.SelectionPointLabelFontFile = ''
clip_partconDisplay.PolarAxes = 'PolarAxesRepresentation'
clip_partconDisplay.ScalarOpacityUnitDistance = 0.11966160496371617

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip_partconDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip_partconDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip_partconDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip_partconDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clip_partconDisplay.DataAxesGrid.XTitleFontFile = ''
clip_partconDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clip_partconDisplay.DataAxesGrid.YTitleFontFile = ''
clip_partconDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clip_partconDisplay.DataAxesGrid.ZTitleFontFile = ''
clip_partconDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
clip_partconDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clip_partconDisplay.DataAxesGrid.XLabelFontFile = ''
clip_partconDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clip_partconDisplay.DataAxesGrid.YLabelFontFile = ''
clip_partconDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clip_partconDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip_partconDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clip_partconDisplay.PolarAxes.PolarAxisTitleFontFile = ''
clip_partconDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clip_partconDisplay.PolarAxes.PolarAxisLabelFontFile = ''
clip_partconDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clip_partconDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
clip_partconDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clip_partconDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from clip_line
clip_lineDisplay = Show(clip_line, renderView1)

# trace defaults for the display properties.
clip_lineDisplay.Representation = 'Surface'
clip_lineDisplay.AmbientColor = [1.0, 0.0, 0.0]
clip_lineDisplay.ColorArrayName = ['POINTS', '']
clip_lineDisplay.DiffuseColor = [0.0, 0.0, 0.0]
clip_lineDisplay.LineWidth = 3.0
clip_lineDisplay.RenderLinesAsTubes = 1
clip_lineDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clip_lineDisplay.SelectOrientationVectors = 'None'
clip_lineDisplay.ScaleFactor = 0.017136803269386294
clip_lineDisplay.SelectScaleArray = 'None'
clip_lineDisplay.GlyphType = 'Arrow'
clip_lineDisplay.GlyphTableIndexArray = 'None'
clip_lineDisplay.GaussianRadius = 0.0008568401634693146
clip_lineDisplay.SetScaleArray = [None, '']
clip_lineDisplay.ScaleTransferFunction = 'PiecewiseFunction'
clip_lineDisplay.OpacityArray = [None, '']
clip_lineDisplay.OpacityTransferFunction = 'PiecewiseFunction'
clip_lineDisplay.DataAxesGrid = 'GridAxesRepresentation'
clip_lineDisplay.SelectionCellLabelFontFile = ''
clip_lineDisplay.SelectionPointLabelFontFile = ''
clip_lineDisplay.PolarAxes = 'PolarAxesRepresentation'
clip_lineDisplay.ScalarOpacityUnitDistance = 0.0368140035535702

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip_lineDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip_lineDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip_lineDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip_lineDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clip_lineDisplay.DataAxesGrid.XTitleFontFile = ''
clip_lineDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clip_lineDisplay.DataAxesGrid.YTitleFontFile = ''
clip_lineDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clip_lineDisplay.DataAxesGrid.ZTitleFontFile = ''
clip_lineDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
clip_lineDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clip_lineDisplay.DataAxesGrid.XLabelFontFile = ''
clip_lineDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clip_lineDisplay.DataAxesGrid.YLabelFontFile = ''
clip_lineDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clip_lineDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip_lineDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clip_lineDisplay.PolarAxes.PolarAxisTitleFontFile = ''
clip_lineDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clip_lineDisplay.PolarAxes.PolarAxisLabelFontFile = ''
clip_lineDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clip_lineDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
clip_lineDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clip_lineDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from clippart
clippartDisplay = Show(clippart, renderView1)

# trace defaults for the display properties.
clippartDisplay.Representation = 'Surface'
clippartDisplay.AmbientColor = [0.0, 0.0, 0.0]
clippartDisplay.ColorArrayName = ['POINTS', '']
clippartDisplay.DiffuseColor = [1.0, 0.4980392156862745, 0.054901960784313725]
clippartDisplay.PointSize = 10.0
clippartDisplay.RenderPointsAsSpheres = 1
clippartDisplay.OSPRayScaleArray = 'sp'
clippartDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clippartDisplay.SelectOrientationVectors = 'None'
clippartDisplay.ScaleFactor = 0.01675102746781083
clippartDisplay.SelectScaleArray = 'sp'
clippartDisplay.GlyphType = 'Arrow'
clippartDisplay.GlyphTableIndexArray = 'sp'
clippartDisplay.GaussianRadius = 0.0008375513733905415
clippartDisplay.SetScaleArray = ['POINTS', 'sp']
clippartDisplay.ScaleTransferFunction = 'PiecewiseFunction'
clippartDisplay.OpacityArray = ['POINTS', 'sp']
clippartDisplay.OpacityTransferFunction = 'PiecewiseFunction'
clippartDisplay.DataAxesGrid = 'GridAxesRepresentation'
clippartDisplay.SelectionCellLabelFontFile = ''
clippartDisplay.SelectionPointLabelFontFile = ''
clippartDisplay.PolarAxes = 'PolarAxesRepresentation'
clippartDisplay.ScalarOpacityUnitDistance = 0.03724402241281036

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clippartDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clippartDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clippartDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clippartDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clippartDisplay.DataAxesGrid.XTitleFontFile = ''
clippartDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clippartDisplay.DataAxesGrid.YTitleFontFile = ''
clippartDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clippartDisplay.DataAxesGrid.ZTitleFontFile = ''
clippartDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
clippartDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clippartDisplay.DataAxesGrid.XLabelFontFile = ''
clippartDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clippartDisplay.DataAxesGrid.YLabelFontFile = ''
clippartDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clippartDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clippartDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clippartDisplay.PolarAxes.PolarAxisTitleFontFile = ''
clippartDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clippartDisplay.PolarAxes.PolarAxisLabelFontFile = ''
clippartDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clippartDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
clippartDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clippartDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from slice
sliceDisplay = Show(slice, renderView1)

# trace defaults for the display properties.
sliceDisplay.Representation = 'Surface'
sliceDisplay.AmbientColor = [0.0, 0.0, 0.0]
sliceDisplay.ColorArrayName = [None, '']
sliceDisplay.Opacity = 0.5
sliceDisplay.RenderLinesAsTubes = 1
sliceDisplay.EdgeColor = [0.0, 0.0, 0.0]
sliceDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
sliceDisplay.SelectOrientationVectors = 'None'
sliceDisplay.ScaleFactor = 0.01675102746781083
sliceDisplay.SelectScaleArray = 'None'
sliceDisplay.GlyphType = 'Arrow'
sliceDisplay.GlyphTableIndexArray = 'None'
sliceDisplay.GaussianRadius = 0.0008375513733905415
sliceDisplay.SetScaleArray = [None, '']
sliceDisplay.ScaleTransferFunction = 'PiecewiseFunction'
sliceDisplay.OpacityArray = [None, '']
sliceDisplay.OpacityTransferFunction = 'PiecewiseFunction'
sliceDisplay.DataAxesGrid = 'GridAxesRepresentation'
sliceDisplay.SelectionCellLabelFontFile = ''
sliceDisplay.SelectionPointLabelFontFile = ''
sliceDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
sliceDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sliceDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sliceDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
sliceDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
sliceDisplay.DataAxesGrid.XTitleFontFile = ''
sliceDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
sliceDisplay.DataAxesGrid.YTitleFontFile = ''
sliceDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
sliceDisplay.DataAxesGrid.ZTitleFontFile = ''
sliceDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
sliceDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
sliceDisplay.DataAxesGrid.XLabelFontFile = ''
sliceDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
sliceDisplay.DataAxesGrid.YLabelFontFile = ''
sliceDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
sliceDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
sliceDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
sliceDisplay.PolarAxes.PolarAxisTitleFontFile = ''
sliceDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
sliceDisplay.PolarAxes.PolarAxisLabelFontFile = ''
sliceDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
sliceDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
sliceDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
sliceDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from clip1
clip1Display = Show(clip1, renderView1)

# trace defaults for the display properties.
clip1Display.Representation = 'Surface With Edges'
clip1Display.AmbientColor = [0.0, 0.0, 0.0]
clip1Display.ColorArrayName = ['POINTS', '']
clip1Display.DiffuseColor = [0.7450980392156863, 0.7450980392156863, 0.7450980392156863]
clip1Display.LineWidth = 0.5
clip1Display.RenderLinesAsTubes = 1
clip1Display.EdgeColor = [0.0, 0.0, 0.0]
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.01720564663410187
clip1Display.SelectScaleArray = 'c'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'c'
clip1Display.GaussianRadius = 0.0008602823317050934
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.SelectionCellLabelFontFile = ''
clip1Display.SelectionPointLabelFontFile = ''
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.10254083672537753

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.XTitleFontFile = ''
clip1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.YTitleFontFile = ''
clip1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.ZTitleFontFile = ''
clip1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.XLabelFontFile = ''
clip1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.YLabelFontFile = ''
clip1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.PolarAxisTitleFontFile = ''
clip1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.PolarAxisLabelFontFile = ''
clip1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from slice2
slice2Display = Show(slice2, renderView1)

# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.AmbientColor = [0.0, 0.0, 0.0]
slice2Display.ColorArrayName = [None, '']
slice2Display.Opacity = 0.5
slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice2Display.SelectOrientationVectors = 'None'
slice2Display.ScaleFactor = 0.035743119076705875
slice2Display.SelectScaleArray = 'None'
slice2Display.GlyphType = 'Arrow'
slice2Display.GlyphTableIndexArray = 'None'
slice2Display.GaussianRadius = 0.0017871559538352936
slice2Display.SetScaleArray = [None, '']
slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
slice2Display.OpacityArray = [None, '']
slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
slice2Display.DataAxesGrid = 'GridAxesRepresentation'
slice2Display.SelectionCellLabelFontFile = ''
slice2Display.SelectionPointLabelFontFile = ''
slice2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice2Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice2Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice2Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
slice2Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
slice2Display.DataAxesGrid.XTitleFontFile = ''
slice2Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
slice2Display.DataAxesGrid.YTitleFontFile = ''
slice2Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
slice2Display.DataAxesGrid.ZTitleFontFile = ''
slice2Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
slice2Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
slice2Display.DataAxesGrid.XLabelFontFile = ''
slice2Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
slice2Display.DataAxesGrid.YLabelFontFile = ''
slice2Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
slice2Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
slice2Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
slice2Display.PolarAxes.PolarAxisTitleFontFile = ''
slice2Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
slice2Display.PolarAxes.PolarAxisLabelFontFile = ''
slice2Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
slice2Display.PolarAxes.LastRadialAxisTextFontFile = ''
slice2Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
slice2Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from extractEdges1
extractEdges1Display = Show(extractEdges1, renderView1)

# trace defaults for the display properties.
extractEdges1Display.Representation = 'Surface'
extractEdges1Display.AmbientColor = [0.0, 0.0, 0.0]
extractEdges1Display.ColorArrayName = [None, '']
extractEdges1Display.DiffuseColor = [0.0, 0.0, 0.0]
extractEdges1Display.RenderLinesAsTubes = 1
extractEdges1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractEdges1Display.SelectOrientationVectors = 'None'
extractEdges1Display.ScaleFactor = 0.03548373579978943
extractEdges1Display.SelectScaleArray = 'None'
extractEdges1Display.GlyphType = 'Arrow'
extractEdges1Display.GlyphTableIndexArray = 'None'
extractEdges1Display.GaussianRadius = 0.0017741867899894715
extractEdges1Display.SetScaleArray = [None, '']
extractEdges1Display.ScaleTransferFunction = 'PiecewiseFunction'
extractEdges1Display.OpacityArray = [None, '']
extractEdges1Display.OpacityTransferFunction = 'PiecewiseFunction'
extractEdges1Display.DataAxesGrid = 'GridAxesRepresentation'
extractEdges1Display.SelectionCellLabelFontFile = ''
extractEdges1Display.SelectionPointLabelFontFile = ''
extractEdges1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
extractEdges1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extractEdges1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extractEdges1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
extractEdges1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
extractEdges1Display.DataAxesGrid.XTitleFontFile = ''
extractEdges1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
extractEdges1Display.DataAxesGrid.YTitleFontFile = ''
extractEdges1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
extractEdges1Display.DataAxesGrid.ZTitleFontFile = ''
extractEdges1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
extractEdges1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
extractEdges1Display.DataAxesGrid.XLabelFontFile = ''
extractEdges1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
extractEdges1Display.DataAxesGrid.YLabelFontFile = ''
extractEdges1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
extractEdges1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
extractEdges1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
extractEdges1Display.PolarAxes.PolarAxisTitleFontFile = ''
extractEdges1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
extractEdges1Display.PolarAxes.PolarAxisLabelFontFile = ''
extractEdges1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
extractEdges1Display.PolarAxes.LastRadialAxisTextFontFile = ''
extractEdges1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
extractEdges1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from extractEdges2
extractEdges2Display = Show(extractEdges2, renderView1)

# trace defaults for the display properties.
extractEdges2Display.Representation = 'Surface'
extractEdges2Display.AmbientColor = [0.0, 0.0, 0.0]
extractEdges2Display.ColorArrayName = [None, '']
extractEdges2Display.DiffuseColor = [0.0, 0.0, 0.0]
extractEdges2Display.RenderLinesAsTubes = 1
extractEdges2Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractEdges2Display.SelectOrientationVectors = 'None'
extractEdges2Display.ScaleFactor = 0.031535136699676516
extractEdges2Display.SelectScaleArray = 'None'
extractEdges2Display.GlyphType = 'Arrow'
extractEdges2Display.GlyphTableIndexArray = 'None'
extractEdges2Display.GaussianRadius = 0.0015767568349838257
extractEdges2Display.SetScaleArray = [None, '']
extractEdges2Display.ScaleTransferFunction = 'PiecewiseFunction'
extractEdges2Display.OpacityArray = [None, '']
extractEdges2Display.OpacityTransferFunction = 'PiecewiseFunction'
extractEdges2Display.DataAxesGrid = 'GridAxesRepresentation'
extractEdges2Display.SelectionCellLabelFontFile = ''
extractEdges2Display.SelectionPointLabelFontFile = ''
extractEdges2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
extractEdges2Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extractEdges2Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extractEdges2Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
extractEdges2Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
extractEdges2Display.DataAxesGrid.XTitleFontFile = ''
extractEdges2Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
extractEdges2Display.DataAxesGrid.YTitleFontFile = ''
extractEdges2Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
extractEdges2Display.DataAxesGrid.ZTitleFontFile = ''
extractEdges2Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
extractEdges2Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
extractEdges2Display.DataAxesGrid.XLabelFontFile = ''
extractEdges2Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
extractEdges2Display.DataAxesGrid.YLabelFontFile = ''
extractEdges2Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
extractEdges2Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
extractEdges2Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
extractEdges2Display.PolarAxes.PolarAxisTitleFontFile = ''
extractEdges2Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
extractEdges2Display.PolarAxes.PolarAxisLabelFontFile = ''
extractEdges2Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
extractEdges2Display.PolarAxes.LastRadialAxisTextFontFile = ''
extractEdges2Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
extractEdges2Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from line
lineDisplay = Show(line, renderView1)

# trace defaults for the display properties.
lineDisplay.Representation = 'Surface'
lineDisplay.ColorArrayName = [None, '']
lineDisplay.DiffuseColor = [0.0, 0.0, 0.0]
lineDisplay.LineWidth = 2.0
lineDisplay.RenderLinesAsTubes = 1
lineDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
lineDisplay.SelectOrientationVectors = 'None'
lineDisplay.ScaleFactor = 0.035743119076705875
lineDisplay.SelectScaleArray = 'None'
lineDisplay.GlyphType = 'Arrow'
lineDisplay.GlyphTableIndexArray = 'None'
lineDisplay.GaussianRadius = 0.0017871559538352936
lineDisplay.SetScaleArray = [None, '']
lineDisplay.ScaleTransferFunction = 'PiecewiseFunction'
lineDisplay.OpacityArray = [None, '']
lineDisplay.OpacityTransferFunction = 'PiecewiseFunction'
lineDisplay.DataAxesGrid = 'GridAxesRepresentation'
lineDisplay.SelectionCellLabelFontFile = ''
lineDisplay.SelectionPointLabelFontFile = ''
lineDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
lineDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
lineDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
lineDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
lineDisplay.DataAxesGrid.XTitleFontFile = ''
lineDisplay.DataAxesGrid.YTitleFontFile = ''
lineDisplay.DataAxesGrid.ZTitleFontFile = ''
lineDisplay.DataAxesGrid.XLabelFontFile = ''
lineDisplay.DataAxesGrid.YLabelFontFile = ''
lineDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
lineDisplay.PolarAxes.PolarAxisTitleFontFile = ''
lineDisplay.PolarAxes.PolarAxisLabelFontFile = ''
lineDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
lineDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(line)
# ----------------------------------------------------------------