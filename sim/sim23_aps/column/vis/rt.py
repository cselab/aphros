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
mf = "m.json"
open(mf, 'w').write('''
{
  "family" : "OSPRay",
  "version" : "0.0",
  "materials" : {
    "water" : {
      "type": "Glass",
      "doubles" : {
          "attenuationColor" : [0.22, 0.34, 0.47],
          "etaInside" : [1.33]
      }
    }
  }
}
''')
materialLibrary1 = GetMaterialLibrary()
materialLibrary1.LoadMaterials = mf

ext = 1

def S(v):
  return np.array(v) * ext


C1 = [
[0.015, 0.027013096538123955, 0.059872780080358526],
[0.015, 0.00709686779203246, 0.008504584363965673],
[0.0, 0.9463305128364528, -0.32320049578349425]
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
renderView1.ViewSize = [1920,1080]
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
ospray = 1
if hasattr(renderView1, 'EnableOSPray'):
    renderView1.EnableOSPRay = ospray
    renderView1.OSPRayRenderer = 'pathtracer'
if hasattr(renderView1, 'EnableRayTracing'):
    renderView1.EnableRayTracing = ospray
    renderView1.BackEnd = 'pathtracer'
renderView1.UseLight = 0
renderView1.AdditionalLights = light1
renderView1.AmbientSamples = 0
renderView1.SamplesPerPixel = 5 if draft else 50
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

surf = Calculator(Input=surf)
surf.ResultNormals = 1
surf.AttributeType = 'Point Data'
surf.ResultArrayName = 'normals'
surf.Function = 'nn'

surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = 'Surface'
surfDisplay.OSPRayMaterial = 'water'


# create a new 'Plane'
planecells = Plane()
eps = 1e-6
planecells.Origin = [-0.02, eps, -0.02]
planecells.Point1 = [0.04, eps, -0.02]
planecells.Point2 = [-0.02, eps, 0.04]
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
