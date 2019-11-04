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
    sys.stderr.write('''usage: {:} [f_*.xmf]
Plots bubbles.
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

# vf input
ff = natsorted(av[1:])
if not len(ff):
    Error("empty file list")
# vf basename
ffb = list(map(os.path.basename, ff))
# vf dirname
ffd = list(map(os.path.dirname, ff))
# steps
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]
# velocity input
ffux = [os.path.join(d, "u.x_{:04d}.xmf".format(s)) for d,s in zip(ffd,ss)]
ffuy = [os.path.join(d, "u.y_{:04d}.xmf".format(s)) for d,s in zip(ffd,ss)]
ffuz = [os.path.join(d, "u.z_{:04d}.xmf".format(s)) for d,s in zip(ffd,ss)]

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

# standard view from corner
C1 = [
 [6.43080086005585, -15.166333096764594, 9.914876482661194],
 [3.0001989267510476, 4.2895772779323025, 2.724262644294233],
 [-0.05939117461388465, 0.33682408883346476, 0.9396926207859086],
 4.2537082476917885,
    ]

C2 = [
  [-10, -10, 3.141592],
  [0, 0, 3.141592],
  [0, 0, 1],
  3.12
    ]

# view in direction (1,1,0)

CC = [C1, C2]
C = CC[cam - 1]


# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2100,1500] if cam == 2 else [1800, 2000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [3.1415927410125732, 3.1415927410125732, 3.1415927410125732]
renderView1.StereoType = 0
renderView1.CameraPosition = C[0]
renderView1.CameraFocalPoint = C[1]
renderView1.CameraViewUp = C[2]
renderView1.CameraParallelScale = C[3]
renderView1.CameraParallelProjection = 1
renderView1.Background = [1., 1., 1.]
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

ext = os.path.splitext(ff[0])[1]
if ext == ".vtk":
    vf = LegacyVTKReader(FileNames=ff)
elif ext == ".xmf":
    vf = XDMFReader(FileNames=ff)
    vf.CellArrayStatus = ['f']
    vf.GridStatus = ['Grid_0']
    ux = XDMFReader(FileNames=ffux)
    ux.CellArrayStatus = ['u.x']
    ux.GridStatus = ['Grid_0']
    uy = XDMFReader(FileNames=ffuy)
    uy.CellArrayStatus = ['u.y']
    uy.GridStatus = ['Grid_0']
    uz = XDMFReader(FileNames=ffuz)
    uz.CellArrayStatus = ['u.z']
    uz.GridStatus = ['Grid_0']
else:
    Error("unknown extension '{:}'".format(ext))

def FT(s):
  vs.append(s)
  vt.append(np.array(s.TimestepValues))
  s = ForceTime(s)
  vft.append(s)
  return s

# sources
vs = []
# time steps
vt = []
# force time
vft = []

vf = FT(vf)
ux = FT(ux)
uy = FT(uy)
uz = FT(uz)

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

append = AppendAttributes(Input=[ux, uy, uz])

cellDatatoPointData1 = CellDatatoPointData(Input=vf)

# create a new 'Contour'
contour1 = Contour(Input=cellDatatoPointData1)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Calculator'
calculator1 = Calculator(Input=append)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'vel'
calculator1.Function = 'iHat*u.x+jHat*u.y+kHat*u.z'

# create a new 'Gradient Of Unstructured DataSet'
grad = GradientOfUnstructuredDataSet(Input=calculator1)
grad.ScalarArray = ['CELLS', 'vel']
grad.ComputeGradient = 0
grad.ComputeVorticity = 1
grad.VorticityArrayName = 'om'


# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# outline
'''
vfDisplay = Show(vf, renderView1)
vfDisplay.Representation = 'Outline'
vfDisplay.ColorArrayName = ['CELLS', '']
vfDisplay.LineWidth = 3.0
vfDisplay.AmbientColor = [0.0, 0.0, 0.0]
'''

# show data from grad
gradDisplay = Show(grad, renderView1)

# get color transfer function/color map for 'omm'
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [1.0, 0.231373, 0.298039, 0.752941, 3.0000000000000004, 0.865003, 0.865003, 0.865003, 5.0, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [1.0, 0.0, 0.5, 0.0, 3.1808218955993652, 0.07236842066049576, 0.5, 0.0, 5.0, 0.2565789520740509, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
gradDisplay.Representation = 'Volume'
gradDisplay.ColorArrayName = ['CELLS', 'om']
gradDisplay.LookupTable = ommLUT
gradDisplay.ScalarOpacityUnitDistance = 0.04251092259923938
gradDisplay.ScalarOpacityFunction = ommPWF

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']


#####################################################
### END OF STATE FILE
#####################################################

SetTime(1)
for i in list(range(len(ss))):
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue

    SetTime(i)

    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1)

exit(0)
