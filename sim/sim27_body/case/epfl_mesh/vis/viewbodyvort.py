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

# vf input
ff = natsorted(av[1:])
# vf basename
ffb = list(map(os.path.basename, ff))
# vf dirname
ffd = list(map(os.path.dirname, ff))
# steps
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]
# omm input
ffomm = [os.path.join(d, "omm_{:04d}.xmf".format(s)) for d,s in zip(ffd,ss)]

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
materialLibrary1 = GetMaterialLibrary()

# view
C1 = [
[2.096370954506405, -3.6052774168001798, 2.482516426770796],
[0.9631715957567765, 0.6238801651579048, 0.44086010664963254],
[-0.10938165494661517, 0.40821789367673533, 0.9063077870366498],
    ]

C2 = [
[2.5308382285072826, 2.109011742127933, 4.655537288725922],
[1.0099864147510529, 0.49055504669887895, 0.477031272094079],
[-0.11697777844051133, 0.9396926207859085, -0.3213938048432697],
    ]

PPS = [0.7829367490651318, 0.753912788055062,]

CC = [C1]
C = CC[cam - 1]

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.OrientationAxesVisibility = 0
renderView1.CameraPosition = C[0]
renderView1.CameraFocalPoint = C[1]
renderView1.CameraViewUp = C[2]
renderView1.Background = [1.0]*3
renderView1.OSPRayMaterialLibrary = materialLibrary1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.CameraParallelScale = PPS[cam - 1]
renderView1.CameraParallelProjection = 1


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
surf = LegacyVTKReader(FileNames=ff)

omm = XDMFReader(FileNames=ffomm)
omm.CellArrayStatus = ['omm']
omm.GridStatus = ['Grid_1']

# list of all sources
vs = [surf, omm]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
surf = ForceTime(surf)
omm = ForceTime(omm)

# all ForceTime
vft = [surf, omm]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

surf = Calculator(Input=surf)
surf.ResultNormals = 1
surf.AttributeType = 'Point Data'
surf.ResultArrayName = 'normals'
surf.Function = 'nn'

def rgb(r, g, b):
    m = 255.
    return [r/m, g/m, b/m]

eps=1e-3
bubbles = Clip(Input=surf)
bubbles.ClipType = 'Box'
bubbles.ClipType.Position = [eps, eps, eps]
bubbles.ClipType.Length = [2 - eps, 1 - eps, 1 - eps]
bubblesDisplay = Show(bubbles, renderView1)
bubblesDisplay.Representation = 'Surface'
bubblesDisplay.ColorArrayName = [None, '']
bubblesDisplay.DiffuseColor = rgb(255, 127, 14)
bubblesDisplay.Ambient = 0.2

bcvtk = LegacyVTKReader(FileNames=['../bc.vtk'])
clip1 = Clip(Input=bcvtk)
clip1.ClipType = 'Box'
clip1.ClipType.Position = [eps, eps, eps]
clip1.ClipType.Length = [2 - eps, 1 - eps, 1 - eps]
clip1Display = Show(clip1, renderView1)
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['POINTS', '']
clip1Display.DiffuseColor = [1.0, 1.0, 1.0]

om0=1
k=0.85
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [om0, 0.231373, 0.298039, 0.752941, 25.0*k, 0.865003, 0.865003, 0.865003, 50.0*k, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [om0, 0.0, 0.5, 0.0, 50.0*k, 1.0, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1

ommvfDisplay = Show(omm, renderView1)
ommvfDisplay.Representation = 'Volume'
ommvfDisplay.AmbientColor = [0.0, 0.0, 0.0]
ommvfDisplay.ColorArrayName = ['CELLS', 'omm']
ommvfDisplay.LookupTable = ommLUT
ommvfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ommvfDisplay.SelectOrientationVectors = 'None'
ommvfDisplay.OpacityArray = [None, '']
ommvfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ommvfDisplay.DataAxesGrid = 'GridAxesRepresentation'
ommvfDisplay.PolarAxes = 'PolarAxesRepresentation'
ommvfDisplay.ScalarOpacityUnitDistance = 0.1
ommvfDisplay.ScalarOpacityFunction = ommPWF

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
