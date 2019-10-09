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


cam = 1
if CheckFlag('-C1'):
    cam = 1
if CheckFlag('-C2'):
    cam = 2
if CheckFlag('-C3'):
    cam = 3
if CheckFlag('-C4'):
    cam = 4

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

light2 = CreateLight()
light2.Intensity = 40.0
light2.Type = 'Positional'
light2.Position = [0.75, 5.0, 5.0]
light2.FocalPoint = [0.75, 2.0, 0.0]
light2.ConeAngle = 15.0
light2.Radius = 0.5


light3 = CreateLight()
light3.Intensity = 30.0
light3.Type = 'Positional'
light3.Position = [0.7500000000000003, 5.0, 2.0]
light3.FocalPoint = [0.7500000000000003, 0.5, 0.75]
light3.ConeAngle = 14.0
light3.Radius = 0.2

if cam == 4:
    light3.Intensity = 25.0
    light3.Type = 'Positional'
    light3.ConeAngle = 15.0
    light3.Radius = 0.3


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


# front
C1 = [
[0.7500000000000001, 0.7516697075897895, 1.9480740875777498],
[0.7500000000000001, 0.6365911653995612, 0.4524949592373693],
[0.0, 0.9970527522269202, -0.07671902812681881]
    ]

# top
C2 = [
[0.7500000000000003, 2.887798073474247, 1.5346513044321786],
[0.7500000000000003, 0.4997215903702232, 0.7860862126109397],
[0.0, 0.3085046033011451, -0.951222849674041]
    ]

# below
C3 = [
[0.7500000000000002, 0.05, 2.7],
[0.7500000000000002, 0.65, 0.4],
[0.0, 0.9825952104266916, -0.1857596631309544]
    ]

# topsingle
C4 = [
    [0.7500000000000006, 2.963107346039059, 2.2607274645150683],
    [0.7500000000000006, 0.5009592363008704, 0.7179205630248001],
    [0.0, 0.5309797150079144, -0.8473845303344368]
        ]

CC = [C1, C2, C3, C4]
C = CC[cam - 1]

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920,1080]
if draft:
    renderView1.ViewSize[0] = renderView1.ViewSize[0] // 2
    renderView1.ViewSize[1] = renderView1.ViewSize[1] // 2
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 0
renderView1.CameraPosition = S(C[0])
renderView1.CameraFocalPoint = S(C[1])
renderView1.CameraViewUp = C[2]

renderView1.Background = [0.0]*3
ospray = 1
if hasattr(renderView1, 'EnableOSPray'):
    renderView1.EnableOSPRay = ospray
    renderView1.OSPRayRenderer = 'pathtracer'
if hasattr(renderView1, 'EnableRayTracing'):
    renderView1.EnableRayTracing = ospray
    renderView1.BackEnd = 'pathtracer'
    renderView1.Denoise = 1
renderView1.UseLight = 0
renderView1.AdditionalLights = light3 if cam == 4 else [light2, light3]
renderView1.AmbientSamples = 1
renderView1.SamplesPerPixel = 5 if draft else 10
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

surf = Transform(Input=surf)
surf.Transform = 'Transform'
surf.Transform.Scale = [50.0]*3

surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = 'Surface'
surfDisplay.OSPRayMaterial = 'water'


bottom = Plane()
bottom.Origin = [-2.25, 0.001, -2.25]
bottom.Point1 = [3.75, 0.001, -2.25]
bottom.Point2 = [-2.25, 0.001, 3.75]
bottom.XResolution = 30
bottom.YResolution = 30
bottomDisplay = Show(bottom, renderView1)
bottomDisplay.Representation = 'Surface With Edges'
bottomDisplay.AmbientColor = [1.]*3
bottomDisplay.ColorArrayName = [None, '']
bottomDisplay.DiffuseColor = [1.]*3
bottomDisplay.LineWidth = 0.5
bottomDisplay.EdgeColor = [0.75]*3


back = Plane()
back.Origin = [-2.25, 0.001, -1.25]
back.Point1 = [3.75, 0.001, -1.25]
back.Point2 = [-2.25, 6.0, -1.25]
back.XResolution = 30
back.YResolution = 30
backDisplay = Show(back, renderView1)
backDisplay.Representation = 'Surface With Edges'
backDisplay.AmbientColor = [1.]*3
backDisplay.ColorArrayName = [None, '']
backDisplay.DiffuseColor = [1.]*3
backDisplay.LineWidth = bottomDisplay.LineWidth
backDisplay.EdgeColor = [0.75]*3


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
