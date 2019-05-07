#!/usr/bin/env pvbatch

# state file generated using paraview version 5.6.0


import sys
import os


av = sys.argv
if len(av) < 2:
    sys.stderr.write('''usage: {:} dat
dat: folder with s_*.vtk, partit_*.csv
'''.format(av[0]))
    exit(1)

# data folder
dat = os.path.abspath(av[1])

def F(f):
    return os.path.join(dat, f)

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()


renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000, 1000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5253897160291672, 0.5222961753606796, 0.5137387216091156]
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.CameraPosition = [0.7660992319552715, 0.2743045036037693, 0.482576870184316]
renderView1.CameraFocalPoint = [0.38399099615374427, 0.667972613582416, 0.53204396344163]
renderView1.CameraViewUp = [0.6613124040614362, 0.5825078997121839, 0.47259967309276063]

cam = F("cam.py")
if os.path.isfile(cam):
    print("Reading custom camera position")
    print(cam)
    ll = open(cam).readlines()
    renderView1.CameraPosition = eval(ll[0])
    renderView1.CameraFocalPoint = eval(ll[1])
    renderView1.CameraViewUp = eval(ll[2])

renderView1.CameraParallelScale = 0.14256872206941143
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableOSPRay = 1
renderView1.AmbientSamples = 1
renderView1.SamplesPerPixel = 50
renderView1.ProgressivePasses = 1
materialLibrary1 = GetMaterialLibrary()
renderView1.OSPRayMaterialLibrary = materialLibrary1

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
dat_line = LegacyVTKReader(FileNames=[F('sp_0000.vtk')])

# create a new 'Legacy VTK Reader'
dat_inter = LegacyVTKReader(FileNames=[F('s_0000.vtk')])

# create a new 'Legacy VTK Reader'
dat_partcon = LegacyVTKReader(FileNames=[F('o.vtk')])

# create a new 'Programmable Filter'
prog_partcon = ProgrammableFilter(Input=dat_partcon)
prog_partcon.Script = """execfile("prog_part.py") """
prog_partcon.RequestInformationScript = ''
prog_partcon.RequestUpdateExtentScript = ''
prog_partcon.PythonPath = ''

# create a new 'Programmable Filter'
prog_line = ProgrammableFilter(Input=dat_line)
prog_line.Script = """execfile("prog_part.py") """
prog_line.RequestInformationScript = ''
prog_line.RequestUpdateExtentScript = ''
prog_line.PythonPath = ''

# create a new 'Clip'
clip_line = Clip(Input=prog_line)
clip_line.ClipType = 'Scalar'
clip_line.Scalars = ['CELLS', 'sc']
clip_line.Value = 0
clip_line.Invert = 0

# create a new 'CSV Reader'
dat_part = CSVReader(FileName=[F('partit_0000.csv')])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=dat_part)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Programmable Filter'
prog_part = ProgrammableFilter(Input=tableToPoints1)
prog_part.Script = """execfile("prog_part.py")"""
prog_part.RequestInformationScript = ''
prog_part.RequestUpdateExtentScript = ''
prog_part.PythonPath = ''

# create a new 'Clip'
clippart = Clip(Input=prog_part)
clippart.ClipType = 'Scalar'
clippart.Scalars = ['POINTS', 'sp']
clippart.Value = 0
clippart.Invert = 0

# create a new 'Programmable Filter'
slice = ProgrammableFilter(Input=prog_part)
slice.OutputDataSetType = 'vtkPolyData'
slice.Script = """import sys
sys.argv = [0.2, 0.2, 0.07]
execfile("prog_polygon.py")  """
slice.RequestInformationScript = ''
slice.RequestUpdateExtentScript = ''
slice.PythonPath = ''

# create a new 'Clip'
clip_partcon = Clip(Input=prog_partcon)
clip_partcon.ClipType = 'Scalar'
clip_partcon.Scalars = ['CELLS', 'sc']
clip_partcon.Value = 0
clip_partcon.Invert = 0

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

# show data from clip_line
clip_lineDisplay = Show(clip_line, renderView1)

# trace defaults for the display properties.
clip_lineDisplay.Representation = 'Surface'
clip_lineDisplay.AmbientColor = [1.0, 0.0, 0.0]
clip_lineDisplay.ColorArrayName = ['POINTS', '']
clip_lineDisplay.DiffuseColor = [0.0, 0.0, 0.0]
clip_lineDisplay.LineWidth = 3.0
clip_lineDisplay.RenderLinesAsTubes = 1

# show data from dat_inter
dat_interDisplay = Show(dat_inter, renderView1)

# trace defaults for the display properties.
dat_interDisplay.Representation = 'Surface With Edges'
dat_interDisplay.AmbientColor = [0.0, 0.0, 0.0]
dat_interDisplay.ColorArrayName = ['POINTS', '']
dat_interDisplay.DiffuseColor = [0.7, 0.7, 0.7]
dat_interDisplay.LineWidth = 0.5
dat_interDisplay.RenderLinesAsTubes = 1
dat_interDisplay.EdgeColor = [0.0, 0.0, 0.0]


# show data from clippart
clippartDisplay = Show(clippart, renderView1)

# trace defaults for the display properties.
clippartDisplay.Representation = 'Surface'
clippartDisplay.AmbientColor = [0.0, 0.0, 0.0]
clippartDisplay.ColorArrayName = ['POINTS', '']
clippartDisplay.DiffuseColor = [1.0, 0.4980392156862745, 0.054901960784313725]
clippartDisplay.PointSize = 10.0
clippartDisplay.RenderPointsAsSpheres = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------

SaveScreenshot("a.png", CompressionLevel=9)
