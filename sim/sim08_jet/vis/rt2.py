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


ospray = 1
vort = 1

av = sys.argv
if len(av) < 2 or av[1] == '-h':
    sys.stderr.write('''usage: {:} [vf_*.xmf]
Plots bubbles.
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

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

# Create a new 'Light'
light3 = CreateLight()
light3.Radius = 3.0
light3.Intensity = 0.1

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
    "glass2" : {
      "type": "Glass",
      "doubles" : {
          "attenuationColor" : [0.22, 0.34, 0.47]
          , "etaInside" : [1.2]
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

W=4000
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
#renderView1.CameraPosition = [0.5, 1.5, 7.]
#renderView1.CameraFocalPoint = [0.5, 1.5, 0.5]
#renderView1.CameraParallelScale = 2.2437468197126655

#renderView1.CameraPosition = [0.5, 1.5, 7.]
#renderView1.CameraFocalPoint = [0.5, 1., 0.5]
#renderView1.CameraParallelScale = 2.2437468197126655
#renderView1.CameraViewUp = [0, 0.86, -0.5]

renderView1.CameraPosition = [0.5, 5.618486180905286, 5.408220702143481]
renderView1.CameraFocalPoint = [0.5, 1.5, 0.5]
renderView1.CameraViewUp = [0.0, 0.766044443118978, -0.6427876096865394]

renderView1.CameraPosition = [2.559243090452642, 3.6914007249843075, 6.157723892186589]
renderView1.CameraFocalPoint = [0.5, 1.5, 0.5]
renderView1.CameraViewUp = [-0.1169777784405109, 0.9396926207859084, -0.3213938048432698]


renderView1.Background = [1.]*3
renderView1.EnableOSPRay = 1
renderView1.OSPRayRenderer = 'pathtracer'
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 200
renderView1.LightScale = 0.8
renderView1.AdditionalLights = light3
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

# create a new 'Box'
box1 = Box()
box1.XLength = 1.
box1.YLength = 3.1
box1.ZLength = 0.05
box1.Center = [0.5, 1.5, -0.025]

# create a new 'Box'
box2 = Box()
box2.XLength = 0.05
box2.YLength = 3.1
box2.ZLength = 1.05
box2.Center = [-0.025, 1.5, 0.475]

# create a new 'Box'
box3 = Box()
box3.XLength = 0.05
box3.YLength = 3.1
box3.ZLength = 1.05
box3.Center = [1.025, 1.5, 0.475]


# create a new 'Box'
box4 = Box()
box4.XLength = 1.
box4.YLength = 0.05
box4.ZLength = 1.
box4.Center = [0.5, -0.025, 0.5]

# create a new 'Box'
box5 = Box()
box5.XLength = 1.
box5.YLength = 0.05
box5.ZLength = 1.
box5.Center = [0.5, 3.025, 0.5]


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

# create a new 'Contour'
contour1 = Contour(Input=clpnt)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

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

# show data from box2
bcolor = [0.62, 0.64, 0.69]

bmat = "None"

# back
box1Display = Show(box1, renderView1)
box1Display.Representation = 'Surface'
box1Display.OSPRayMaterial = bmat
box1Display.DiffuseColor = bcolor

# left
box2Display = Show(box2, renderView1)
box2Display.Representation = 'Surface'
box2Display.OSPRayMaterial = bmat
box2Display.DiffuseColor = bcolor

# right
box3Display = Show(box3, renderView1)
box3Display.Representation = 'Surface'
box3Display.OSPRayMaterial = 'glass2'
box3Display.DiffuseColor = bcolor

# bottom
box4Display = Show(box4, renderView1)
box4Display.Representation = 'Surface'
box4Display.OSPRayMaterial = bmat
box4Display.DiffuseColor = bcolor

#box5Display = Show(box5, renderView1)
#box5Display.Representation = 'Surface'
#box5Display.OSPRayMaterial = bmat
#box5Display.DiffuseColor = bcolor2

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
    q = 2 ** 0.5
    q = 2
    w = 4000 // 4
    #SaveScreenshot(fn, renderView1, ImageResolution=[w, int(q * w)])
    SaveScreenshot(fn, renderView1)

exit(0)
