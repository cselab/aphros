#!/usr/bin/env pvbatch

# Le Corbusier colors
# https://www.lescouleurs.ch/fileadmin/media/le-corbusier/polychromie-architecturale/le-corbusier-polychromie-architecturale-1959-1.png

# pass 'draft' as first argument for lowres

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
    sys.stderr.write('''usage: {:} [vf_*.xmf]
Plots bubbles.
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

av = av[1:]

draft = 0
if av[0] == "draft":
  draft = 1
  av = av[1:]
vol = 1
cont = 0

# vf input
ff = natsorted(av[0:])
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
    ,
    "water2" : {
      "type": "Glass",
      "doubles" : {
          "attenuationColor" : [0.22, 0.34, 0.47],
          "attenuationDistance" : [5.],
          "etaInside" : [1.15]
      }
    }
    ,
    "glass2" : {
      "type": "Glass",
      "doubles" : {
          "attenuationColor" : [0.22, 0.34, 0.47]
          , "attenuationDistance" : [5.]
          , "etaInside" : [1.25]
      }
    }
    ,
    "glass" : {
      "type": "Glass"
    }
    ,
    "aluminum" : {
      "type" : "Metal"
    }
  }
}
''')

materialLibrary1 = GetMaterialLibrary()
materialLibrary1.LoadMaterials = mf

W = 2000 if draft else 6000
q=2**0.5
H=int(W*q)

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [W, H]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 1.5, 0.5]
renderView1.UseLight = 0
renderView1.StereoType = 0

if 1:
  renderView1.CameraPosition = [2.55, 3.69, 6.15]
  renderView1.CameraFocalPoint = [0.5, 1.5, 0.5]
  renderView1.CameraViewUp = [-0.11, 0.93, -0.32]
  renderView1.Background = [1.]*3
else:
  renderView1.CameraPosition = [0.5, 1.5, 7]
  renderView1.CameraFocalPoint = [0.5, 1.5, 0.5]
  renderView1.CameraViewUp = [0., 1., 0.]
  renderView1.Background = [0.]*3

if 0:
  renderView1.CameraParallelScale = 1.5
  renderView1.CameraParallelProjection = 1

renderView1.EnableOSPRay = 1
renderView1.OSPRayRenderer = 'pathtracer'
renderView1.AmbientSamples = 0
renderView1.SamplesPerPixel = 200 if draft else 400
renderView1.LightScale = 1.

D = 10.
RAD = 20.
LI = 0.8
light1 = CreateLight()
light1.Radius = RAD
light1.Intensity = LI
light1.FocalPoint = [0]*3
light1.Position = [D,0,0]
#light1.DiffuseColor = [1,0,0]

light2 = CreateLight()
light2.Radius = RAD
light2.Intensity = LI
light2.FocalPoint = [0]*3
light2.Position = [0,D,0]
#light2.DiffuseColor = [0,1,0]

# Create a new 'Light'
light3 = CreateLight()
light3.Radius = RAD
light3.Intensity = LI
light3.FocalPoint = [0]*3
light3.Position = [0,0,D]
#light3.DiffuseColor = [0.,0.,1]


renderView1.AdditionalLights = [light1, light2, light3]
#renderView1.AdditionalLights = [light3]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# BEGIN BOXES
# ----------------------------------------------------------------

eps = -1e-4

# back
box1 = Box()
box1.XLength = 1.
box1.YLength = 3.1
box1.ZLength = 0.05
box1.Center = [0.5, 1.5, -0.025 - eps]

# left
box2 = Box()
box2.XLength = 0.05
box2.YLength = 3.1
box2.ZLength = 1.05
box2.Center = [-0.025 - eps, 1.5, 0.475]

# right
box3 = Box()
box3.XLength = 0.05
box3.YLength = 3.1
box3.ZLength = 1.05
box3.Center = [1.025 + eps, 1.5, 0.475]

# bottom
box4 = Box()
box4.XLength = 1.
box4.YLength = 0.05
box4.ZLength = 1.
box4.Center = [0.5, -0.025 - eps, 0.5]

# top
box5 = Box()
box5.XLength = 1.
box5.YLength = 0.05
box5.ZLength = 1.
box5.Center = [0.5, 3.025 + eps, 0.5]


# ----------------------------------------------------------------
# END BOXES
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
vf = XDMFReader(FileNames=ff)
vf.CellArrayStatus = ['vf']
vf.GridStatus = ['Grid_0']

# list of all sources
vs = [vf]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
vf = ForceTime(vf)

# all ForceTime
vft = [vf]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# create a new 'Cell Data to Point Data'
clpnt = CellDatatoPointData(Input=vf)

rsmp = ResampleToImage(Input=clpnt)
#nx = 128
#nx = 256
nx = 320
rsmp.SamplingDimensions = [nx + 1,  nx * 3 + 1, nx + 1]
rsmp.SamplingBounds = [0.0, 1.0, 0.0, 3.0, 0.0, 1.0]

clpnt = rsmp

if cont:
  # create a new 'Contour'
  contour1 = Contour(Input=clpnt)
  contour1.ContourBy = ['POINTS', 'vf']
  contour1.Isosurfaces = [0.5]
  contour1.PointMergeMethod = 'Uniform Binning'

if vol:
  isoVolume1 = IsoVolume(Input=clpnt)
  isoVolume1.InputScalars = ['POINTS', 'vf']
  isoVolume1.ThresholdRange = [0.0, 0.5]


# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

if cont:
  # show data from contour1
  contour1Display = Show(contour1, renderView1)

  # trace defaults for the display properties.
  contour1Display.Representation = 'Surface'
  contour1Display.ColorArrayName = [None, '']
  contour1Display.LineWidth = 3.0
  contour1Display.OSPRayScaleArray = 'Normals'
  contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
  contour1Display.OSPRayMaterial = 'water'
  contour1Display.SelectOrientationVectors = 'None'
  contour1Display.ScaleFactor = 0.2871810391545296
  contour1Display.SelectScaleArray = 'None'
  contour1Display.GlyphType = 'Arrow'
  contour1Display.GlyphTableIndexArray = 'None'
  contour1Display.GaussianRadius = 0.014359051957726479
  contour1Display.SetScaleArray = ['POINTS', 'Normals']
  contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
  contour1Display.OpacityArray = ['POINTS', 'Normals']
  contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
  contour1Display.DataAxesGrid = 'GridAxesRepresentation'
  contour1Display.SelectionCellLabelFontFile = ''
  contour1Display.SelectionPointLabelFontFile = ''
  contour1Display.PolarAxes = 'PolarAxesRepresentation'

  # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
  contour1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

  # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
  contour1Display.ScaleTransferFunction.Points = [-0.9999940991401672, 1.0, 0.5, 0.0, 0.9999958872795105, 1.0, 0.5, 0.0]

  # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
  contour1Display.OpacityTransferFunction.Points = [-0.9999940991401672, 1.0, 0.5, 0.0, 0.9999958872795105, 1.0, 0.5, 0.0]


if vol:
  # show data from isoVolume1
  isoVolume1Display = Show(isoVolume1, renderView1)

  # trace defaults for the display properties.
  isoVolume1Display.Representation = 'Surface'
  isoVolume1Display.ColorArrayName = ['POINTS', '']
  isoVolume1Display.OSPRayScaleFunction = 'PiecewiseFunction'
  isoVolume1Display.OSPRayMaterial = 'water2'
  isoVolume1Display.SelectOrientationVectors = 'None'
  isoVolume1Display.ScaleFactor = 0.30000000000000004
  isoVolume1Display.SelectScaleArray = 'vf'
  isoVolume1Display.GlyphType = 'Arrow'
  isoVolume1Display.GlyphTableIndexArray = 'vf'
  isoVolume1Display.GaussianRadius = 0.015
  isoVolume1Display.SetScaleArray = [None, '']
  isoVolume1Display.ScaleTransferFunction = 'PiecewiseFunction'
  isoVolume1Display.OpacityArray = [None, '']
  isoVolume1Display.OpacityTransferFunction = 'PiecewiseFunction'
  isoVolume1Display.DataAxesGrid = 'GridAxesRepresentation'
  isoVolume1Display.SelectionCellLabelFontFile = ''
  isoVolume1Display.SelectionPointLabelFontFile = ''
  isoVolume1Display.ScalarOpacityUnitDistance = 0.01796577493112287


# colors
CR = [0.718, 0.204, 0.035]
CG = [0.718, 0.71, 0.0]
CB = [0.192, 0.361, 0.655]
CS = [0.6, 0.6, 0.6]

matn = "None"
matg = "glass2"

boxgl = 1
boxcl = 1

if boxcl:
  # back
  box1Display = Show(box1, renderView1)
  box1Display.Representation = 'Surface'
  box1Display.OSPRayMaterial = matn
  box1Display.DiffuseColor = CB

if boxcl:
  # left
  box2Display = Show(box2, renderView1)
  box2Display.Representation = 'Surface'
  box2Display.OSPRayMaterial = matn
  box2Display.DiffuseColor = CR

if boxgl:
  # right
  box3Display = Show(box3, renderView1)
  box3Display.Representation = 'Surface'
  box3Display.OSPRayMaterial = matg
  box3Display.DiffuseColor = CS

if boxcl:
  # bottom
  box4Display = Show(box4, renderView1)
  box4Display.Representation = 'Surface'
  box4Display.OSPRayMaterial = matn
  box4Display.DiffuseColor = CG

if boxcl:
  pass
  # top
  #box5Display = Show(box5, renderView1)
  #box5Display.Representation = 'Surface'
  #box5Display.OSPRayMaterial = matg
  #box5Display.DiffuseColor = CS

# ----------------------------------------------------------------
# finally, restore active source
#SetActiveSource(box4)
# ----------------------------------------------------------------

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
