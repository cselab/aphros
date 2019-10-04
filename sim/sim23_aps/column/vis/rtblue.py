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
light1.Intensity = 1.35
light1.Position = [0.5, 20.0, 5.0]
light1.FocalPoint = [0.5, 0.0, 0.5]
light1.Radius = 10.0

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
          "attenuationDistance" : [3.0],
          "eta" : [1.33]
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

# top
C2 = [
[0.7500000000000001, 2.8113238646327012, 1.2945544554455446],
[0.7500000000000001, 0.4467885628015819, 0.75],
[0.0, 0.27177120267322913, -0.9623618931553486]
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
renderView1.EnableOSPRay = 1
renderView1.OSPRayRenderer = 'pathtracer'
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
planecells.Origin = [-1.25, 0.001, -1.25]
planecells.Point1 = [2.75, 0.001, -1.25]
planecells.Point2 = [-1.25, 2.75, -1.25]
planecells.XResolution = 20
planecells.YResolution = 20
planecellsDisplay = Show(planecells, renderView1)
planecellsDisplay.Representation = 'Surface With Edges'
planecellsDisplay.AmbientColor = [0.0, 0.0, 0.0]
planecellsDisplay.ColorArrayName = [None, '']
planecellsDisplay.DiffuseColor = [1.]*3
planecellsDisplay.LineWidth = 0.5
planecellsDisplay.EdgeColor = [0.6509803921568628, 0.6509803921568628, 0.6509803921568628]


planevert = Plane()
planevert.Origin = [-1.25, 0.001, -1.25]
planevert.Point1 = [2.75, 0.001, -1.25]
planevert.Point2 = [-1.25, 0.001, 2.75]
planevert.XResolution = planecells.XResolution
planevert.YResolution = planecells.YResolution
planevertDisplay = Show(planevert, renderView1)
planevertDisplay.Representation = 'Surface With Edges'
planevertDisplay.AmbientColor = [0.0, 0.0, 0.0]
planevertDisplay.ColorArrayName = [None, '']
planevertDisplay.DiffuseColor = planecellsDisplay.DiffuseColor
planevertDisplay.LineWidth = planecellsDisplay.LineWidth
planevertDisplay.EdgeColor = planecellsDisplay.EdgeColor


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
