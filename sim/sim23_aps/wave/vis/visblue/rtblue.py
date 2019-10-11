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
if CheckFlag('-C1'):
    cam = 1
if CheckFlag('-C2'):
    cam = 2
if CheckFlag('-C3'):
    cam = 3
if CheckFlag('-C4'):
    cam = 4

draft = CheckFlag('-draft')
fine = CheckFlag('-fine')
fine2 = CheckFlag('-fine2')

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
light2.Intensity = 60.0
light2.Type = 'Positional'
light2.Position = [2.0, 8.0, 1.0]
light2.FocalPoint = [1.0, 1.0, 0.5]
light2.ConeAngle = 11.0
light2.Radius = 0.5

light11 = CreateLight()
light11.Intensity = 100.0
light11.Type = 'Positional'
light11.Position = [2.0, 20.0, 2.0]
light11.FocalPoint = [1.0, 0.0, 0.5]
light11.ConeAngle = 23.0
light11.Radius = 5.0


# get the material library
mf = "m.json"
tp = "checker_grey10x10.png"
open(mf, 'w').write('''
{
  "family" : "OSPRay",
  "version" : "0.0",
  "materials" : {
    "water" : {
      "type": "Glass",
      "doubles" : {
          "attenuationColor" : [0.22, 0.34, 0.47],
          "attenuationDistance" : [1.5],
          "eta" : [1.33]
      }
    },
    "checker grey" : {
      "type" : "OBJMaterial",
      "textures" : {
        "map_d" : "CHECKERGRAY"
      }
    }
  }
}
'''.replace("CHECKERGRAY", os.path.abspath(tp)))
materialLibrary1 = GetMaterialLibrary()
materialLibrary1.LoadMaterials = mf

ext = 1

def S(v):
  return np.array(v) * ext


# view
C1 = [
[1.673947199541654, 1.155534385955274, 3.017460244409912],
[0.49676167169496643, -0.25002651930486053, -1.3449860440897603],
[-0.0893981734829278, 0.9547941520217919, -0.2835067791833273]
    ]

CC = [C1]
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

#renderView1.Background = [0.0]*3
renderView1.Background = [0.2]*3
ospray = 1
if hasattr(renderView1, 'EnableOSPray'):
    renderView1.EnableOSPRay = ospray
    renderView1.OSPRayRenderer = 'pathtracer'
if hasattr(renderView1, 'EnableRayTracing'):
    renderView1.EnableRayTracing = ospray
    renderView1.BackEnd = 'OSPRay pathtracer'
    renderView1.Denoise = 1
renderView1.UseLight = 0
renderView1.AdditionalLights = [light2, light11]
renderView1.LightScale=1
renderView1.AmbientSamples = 0
renderView1.SamplesPerPixel = 5 if draft else (400 if fine2 else (100 if fine else 50))
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

floor = Plane()
floor.Origin = [-20.0, -0.01, -35.0]
floor.Point1 = [-20.0, -0.01, 5.0]
floor.Point2 = [20.0, -0.01, -35.0]
floorDisplay = Show(floor, renderView1)
floorDisplay.Representation = 'Surface'
floorDisplay.ColorArrayName = ['POINTS', '']
floorDisplay.OSPRayMaterial = 'checker grey'

bottom = Plane()
bottom.Origin = [0.0, 0.00015, 0.0]
bottom.Point1 = [2.0, 0.00015, 0.0]
bottom.Point2 = [0.0, 0.00015, 1.0]
bottomDisplay = Show(bottom, renderView1)
bottomDisplay.Representation = 'Surface'
bottomDisplay.ColorArrayName = [None, '']


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
