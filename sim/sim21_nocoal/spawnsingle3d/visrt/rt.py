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
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

def CheckFlag(name):
    if name in av:
        av.remove(name)
        return True
    return False


cam = 1  # view from side perspective
if CheckFlag('-C2'):  # view from top
    cam = 2
if CheckFlag('-C3'):  # view from side parallel
    cam = 3

draft = CheckFlag('-draft')

# vf input
ff = natsorted(av[1:])
# vf basename
ffb = list(map(os.path.basename, ff))
# vf dirname
ffd = list(map(os.path.dirname, ff))
# steps
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]


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

light1 = CreateLight()
light1.Position = [0.01, 0.2, 0.1]
light1.FocalPoint = [0.01, 0.0, 0.01]
light1.Radius = 10
light1.Type = 'Directional'
light1.Intensity = 1.2

# get the material library
materialLibrary1 = GetMaterialLibrary()

ext = 1

def S(v):
  return np.array(v) * ext

C1 = [
[0.01, 0.017161797980427655, 0.05424638773025392],
[0.01, -0.01411619351504851, -0.10666492902928304],
[0.0, 0.981627183447664, -0.19080899537654494],
    ]

C2 = [
[0.01, 0.04969376323491339, 0.02307639700252371],
[0.01, -0.001651995287112607, 0.004388069244686046],
[0.0, 0.3420201433256689, -0.9396926207859084]
    ]

C3 = [
[0.009999999776482582, 0.01, 0.06999441558868025],
[0.009999999776482582, 0.01, 0.009999999776482582],
[0.0, 1.0, 0.0]
    ]

C = C2 if cam == 2 else C3 if cam == 3 else C1

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000,1000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 0
renderView1.CameraPosition = S(C[0])
renderView1.CameraFocalPoint = S(C[1])
renderView1.CameraViewUp = C[2]

if cam == 3:
  renderView1.CameraParallelScale = 0.01*ext
  renderView1.CameraParallelProjection = 1

renderView1.Background = [0.0]*3
renderView1.EnableOSPRay = 1
renderView1.OSPRayRenderer = 'pathtracer'
renderView1.UseLight = 0
renderView1.AdditionalLights = light1
renderView1.AmbientSamples = 0
renderView1.SamplesPerPixel = 5 if draft else 100
renderView1.OSPRayMaterialLibrary = materialLibrary1


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
surf = LegacyVTKReader(FileNames=ff)

# list of all sources
vs = [surf]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
surf = ForceTime(surf)

# all ForceTime
vft = [surf]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

surf = GenerateSurfaceNormals(Input=surf)

clip1 = Clip(Input=surf)
clip1.ClipType = 'Plane'
clip1.Crinkleclip = 1
clip1.ClipType.Origin = [0.009999999776482582, 0.019, 0.009999999776482582]
clip1.ClipType.Normal = [0.0, 1.0, 0.0]
surf=clip1

surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = 'Surface'
surfDisplay.OSPRayMaterial = 'glass'


# create a new 'Plane'
planecells = Plane()
planecells.Origin = [-0.02, 0.0, -0.02]
planecells.Point1 = [0.04, 0.0, -0.02]
planecells.Point2 = [-0.02, 0.0, 0.04]
planecells.XResolution = 15
planecells.YResolution = 15
planecellsDisplay = Show(planecells, renderView1)
planecellsDisplay.Representation = 'Surface With Edges'
planecellsDisplay.AmbientColor = [0.0, 0.0, 0.0]
planecellsDisplay.ColorArrayName = [None, '']
planecellsDisplay.DiffuseColor = [0.8]*3
planecellsDisplay.LineWidth = 0.5
planecellsDisplay.EdgeColor = [0.6509803921568628, 0.6509803921568628, 0.6509803921568628]


planevert = Plane()
planevert.Origin = [-0.02, -0.02, -0.02]
planevert.Point1 = [0.04, -0.02, -0.02]
planevert.Point2 = [-0.02, 0.04, -0.02]
planevert.XResolution = 15
planevert.YResolution = 15
planevertDisplay = Show(planevert, renderView1)
planevertDisplay.Representation = 'Surface With Edges'
planevertDisplay.AmbientColor = [0.0, 0.0, 0.0]
planevertDisplay.ColorArrayName = [None, '']
planevertDisplay.DiffuseColor = planecellsDisplay.DiffuseColor
planevertDisplay.LineWidth = 0.5
planevertDisplay.EdgeColor = [0.6509803921568628, 0.6509803921568628, 0.6509803921568628]


'''
planeglass = Plane()
planeglass.Origin = [-0.02, 0.01488, -0.02]
planeglass.Point1 = [0.04, 0.01488, -0.02]
planeglass.Point2 = [-0.02, 0.01488, 0.04]
planeglass.XResolution = 100
planeglass.YResolution = 100
calculator1 = Calculator(Input=planeglass)
calculator1.AttributeType = 'Point Data'
calculator1.ResultArrayName = 'd'
calculator1.Function = 'max(abs(coordsX-0.01),abs(coordsZ-0.01))'
planeglassclip = Clip(Input=calculator1)
planeglassclip.ClipType = 'Scalar'
planeglassclip.Scalars = ['POINTS', 'd']
planeglassclip.Value = 0.01
planeglassclip.Invert = 0
planeglassclipDisplay = Show(planeglassclip, renderView1)
planeglassclipDisplay.Representation = 'Surface'
planeglassclipDisplay.OSPRayMaterial = 'glass'
'''

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

exit(0)
