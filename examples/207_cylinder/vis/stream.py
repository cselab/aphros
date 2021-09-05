#!/usr/bin/env pvbatch

from paraview.simple import *

from glob import glob
import sys
import re
import math
import os
import numpy as np

def Log(s):
    s += "\n"
    o = sys.stderr
    o.write(s)
    o.flush()

# natural sort
def natkey(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]

def natsorted(v):
    return sorted(v, key=natkey)
# Returns sorted list of files in base by pattern pre_*.xmf

# Sets time of datasets to step i
def SetTime(i):
    global vft, vt
    for j in range(len(vft)):
        s = vft[j]
        s.ForcedTime = vt[j][i]
        s.UpdatePipeline()

# Returns bounding box of object o
def GetBox(o):
    o.UpdatePipeline()
    di = o.GetDataInformation()
    lim = di.DataInformation.GetBounds()
    lim0 = np.array(lim[::2])
    lim1 = np.array(lim[1::2])
    return lim0, lim1


av = sys.argv
if len(av) < 2 or av[1] == '-h':
    sys.stderr.write('''usage: {:} [sm_*.vtk]
Plots isosurface.
# Output:
# a_*.png in current directory
'''.format(av[0]))
    exit(1)

def CheckFlag(name):
    if name in av:
        av.remove(name)
        return True
    return False


cam = 1  # view from side perspective
if CheckFlag('-C1'):
    cam = 1
if CheckFlag('-C2'):
    cam = 2

draft = CheckFlag('-draft')

# vx input
ff = natsorted(av[1:])
# sm basename
ffb = list(map(os.path.basename, ff))
# sm dirname
ffd = list(map(os.path.dirname, ff))
# steps
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]
# omm input
ffvy = [os.path.join(d, "vy_{:04d}.xmf".format(s)) for d,s in zip(ffd,ss)]


# output pattern (:0 substituted by frame number)
bo = "a_{:}.png"

#####################################################
### BEGIN OF STATE FILE
#####################################################

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------

materialLibrary1 = GetMaterialLibrary()


# standard view from corner
C1 = [
    [9.544012387968593, 4.987782847782017, 43.20729827964174],
    [9.544012387968593, 4.987782847782017, 0.009765625],
    [0.0, 1.0, 0.0],
    2.707925050119619,
    ]

CC = [C1]
C = CC[cam - 1]

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2000, 1000]
if draft:
    renderView1.ViewSize[0] = renderView1.ViewSize[0] // 2
    renderView1.ViewSize[1] = renderView1.ViewSize[1] // 2
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 0
renderView1.CameraPosition = C[0]
renderView1.CameraFocalPoint = C[1]
renderView1.CameraViewUp = C[2]
renderView1.CameraParallelScale = C[3]
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.]*3
renderView1.UseLight = 0
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.OSPRayMaterialLibrary = materialLibrary1


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

vx = XDMFReader(FileNames=ff)
vx.CellArrayStatus = ['vx']
vx.GridStatus = ['Grid_1']

vy = XDMFReader(FileNames=ffvy)
vy.CellArrayStatus = ['vy']
vy.GridStatus = ['Grid_1']

# list of all sources
vs = [vx, vy]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
vx = ForceTime(vx)
vy = ForceTime(vy)

# all ForceTime
vft = [vx, vy]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

ebvtk = "../eb.vtk"
if os.path.isfile(ebvtk):
    ebvtk = LegacyVTKReader(FileNames=[ebvtk])
    ebvtk = Clip(Input=ebvtk)
    ebvtk.ClipType = 'Scalar'
    ebvtk.Scalars = ['CELLS', 'face']
    ebvtk.Value = 0.5
    ebvtkDisplay = Show(ebvtk, renderView1)
    ebvtkDisplay.Representation = 'Wireframe'
    ebvtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
    ebvtkDisplay.ColorArrayName = ['POINTS', '']
    ebvtkDisplay.DiffuseColor = [0.0, 0.0, 0.0]
    ebvtkDisplay.LineWidth = 5.

bcvtk = "../bc.vtk"
if os.path.isfile(bcvtk):
    bcvtk = LegacyVTKReader(FileNames=[bcvtk])
    bcvtkDisplay = Show(bcvtk, renderView1)
    bcvtkDisplay.Representation = 'Wireframe'
    bcvtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
    bcvtkDisplay.ColorArrayName = ['POINTS', '']
    bcvtkDisplay.DiffuseColor = [0.0, 0.0, 0.0]
    bcvtkDisplay.LineWidth = 5.


# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(Input=[vx, vy])

# create a new 'Calculator'
calculator1 = Calculator(Input=appendAttributes1)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'vel'
calculator1.Function = 'vx*iHat+vy*jHat+0*kHat'

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=calculator1,
    SeedType='Point Source')
streamTracer1.Vectors = ['CELLS', 'vel']
streamTracer1.SurfaceStreamlines = 1
streamTracer1.MaximumStreamlineLength = 20.0

# init the 'Point Source' selected for 'SeedType'
streamTracer1.SeedType.Center = [7.997848606294397, 5.0301943135265015, 0.009765625]
streamTracer1.SeedType.NumberOfPoints = 1000
streamTracer1.SeedType.Radius = 1.1370834313328775

# show data from streamTracer1
streamTracer1Display = Show(streamTracer1, renderView1)
streamTracer1Display.Representation = 'Wireframe'
streamTracer1Display.AmbientColor = [0.0, 0.0, 0.0]
streamTracer1Display.ColorArrayName = [None, '']

#####################################################
### END OF STATE FILE
#####################################################

for i in list(range(len(ss))):
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue

    SetTime(i)

    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1)
