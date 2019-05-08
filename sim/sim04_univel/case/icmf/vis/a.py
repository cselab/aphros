#!/usr/bin/env pvbatch

# state file generated using paraview version 5.6.0


import sys
import os
import glob


av = sys.argv
if len(av) < 3:
    sys.stderr.write('''usage: {:} dat out
dat: folder with s_*.vtk, partit_*.csv
out: path to output png
'''.format(av[0]))
    exit(1)


# data folder
dat = os.path.abspath(av[1])
# output
out = os.path.abspath(av[2])

# parameters
pinter = True
ppart = True

def F(f):
    return os.path.join(dat, f)

def FF(f):
    return sorted(glob.glob(os.path.join(dat, f)))

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

light2 = CreateLight()
light2.Type = 'Positional'
light2.Position = [0.18546666236455242, -0.6018346614839948, 0.33186906365471086]
light2.FocalPoint = [0.5082919864137087, 0.28176814443494036, 0.48899205815810054]

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000, 1000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5253897160291672, 0.5222961753606796, 0.5137387216091156]
renderView1.UseLight = 0
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.AdditionalLights = light2
renderView1.CameraPosition = [0.24839119124607448, -0.10623063113324308, 0.3704770501303224]
renderView1.CameraFocalPoint = [0.7212417088961819, 0.6928833886588502, 0.5883919766652376]
renderView1.CameraViewUp = [0.7866517756373992, -0.5447244857787051, 0.2906100797970575]

cam = F("cam.py")
if os.path.isfile(cam):
    print("Reading custom camera position")
    print(cam)
    ll = open(cam).readlines()
    renderView1.CameraPosition = eval(ll[0])
    renderView1.CameraFocalPoint = eval(ll[1])
    renderView1.CameraViewUp = eval(ll[2])

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
dat_inter = LegacyVTKReader(FileNames=FF('s_*.vtk'))

# create a new 'Legacy VTK Reader'
dat_partcon = LegacyVTKReader(FileNames=FF('partit_*.vtk'))

# create a new 'Programmable Filter'
prog_partcon = ProgrammableFilter(Input=dat_partcon)
prog_partcon.Script = """execfile("prog_part.py") """
prog_partcon.RequestInformationScript = ''
prog_partcon.RequestUpdateExtentScript = ''
prog_partcon.PythonPath = ''

# create a new 'CSV Reader'
dat_part = CSVReader(FileName=FF('partit_*.csv'))

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
clip1 = Clip(Input=dat_inter)
clip1.ClipType = 'Plane'
clip1.Scalars = ['CELLS', 'c']
clip1.Value = 8007000.0
clip1.Crinkleclip = 1
clip1.ClipType.Origin = [0.5788521482325462, 0.408296197740581, 0.5227404940354464]
clip1.ClipType.Normal = [-0.4573763428286224, 0.8819643456684312, -0.11377949723201745]

# create a new 'Clip'
clippart = Clip(Input=prog_part)
clippart.ClipType = 'Scalar'
clippart.Scalars = ['POINTS', 'sp']
clippart.Value = 0
clippart.Invert = 0

# create a new 'Clip'
clip_partcon = Clip(Input=prog_partcon)
clip_partcon.ClipType = 'Scalar'
clip_partcon.Scalars = ['CELLS', 'sc']
clip_partcon.Value = 0
clip_partcon.Invert = 0

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from clip1
if pinter:
    clip1Display = Show(dat_inter, renderView1)
    clip1Display.Representation = 'Surface With Edges'
    clip1Display.AmbientColor = [0.0, 0.0, 0.0]
    clip1Display.ColorArrayName = ['POINTS', '']
    clip1Display.DiffuseColor = [0.75, 0.75, 0.75]
    clip1Display.LineWidth = 0.5
    clip1Display.RenderLinesAsTubes = 1
    clip1Display.EdgeColor = [0.0, 0.0, 0.0]

if ppart:
    clippartDisplay = Show(clippart, renderView1)
    clippartDisplay.Representation = 'Surface'
    clippartDisplay.AmbientColor = [0.0, 0.0, 0.0]
    clippartDisplay.ColorArrayName = ['POINTS', '']
    clippartDisplay.DiffuseColor = [1.0, 0.4980392156862745, 0.054901960784313725]
    clippartDisplay.PointSize = 10.0
    clippartDisplay.RenderPointsAsSpheres = 1

    clip_partconDisplay = Show(clip_partcon, renderView1)
    clip_partconDisplay.Representation = 'Surface'
    clip_partconDisplay.AmbientColor = [1.0, 0.0, 0.0]
    clip_partconDisplay.ColorArrayName = ['POINTS', '']
    clip_partconDisplay.DiffuseColor = [1.0, 0.4980392156862745, 0.054901960784313725]
    clip_partconDisplay.LineWidth = 2
    clip_partconDisplay.RenderLinesAsTubes = 1

SetActiveSource(None)

sc = GetAnimationScene()
sc.NumberOfFrames = len(dat_inter.TimestepValues)
sc.PlayMode = 'Snap To TimeSteps'
#SaveAnimation(out, FrameWindow=[0,10], CompressionLevel=9)
SaveAnimation(out, CompressionLevel=9)
