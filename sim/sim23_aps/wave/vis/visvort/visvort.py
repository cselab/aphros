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

draft = CheckFlag('-draft')

# sm input
ff = natsorted(av[1:])
# sm basename
ffb = list(map(os.path.basename, ff))
# sm dirname
ffd = list(map(os.path.dirname, ff))
# steps
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]
# vf input
ffvf = [os.path.join(d, "vf_{:04d}.xmf".format(s)) for d,s in zip(ffd,ss)]
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

materialLibrary1 = GetMaterialLibrary()


# view
C1 = [
        [1.0, 0.21919698946269434, 2.0664402974005656],
        [1.0, 1.0409089888243557, -2.5937200255409447],
        [0.0, 0.9848077530122081, 0.17364817766693033]
    ]

# front
C2 = [
        [1.0, 0.5, 3.122550182722647],
        [1.0, 0.5, -1.6095006248462298],
        [0., 1., 0.]
    ]

CC = [C1, C2]
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
renderView1.CameraPosition = C[0]
renderView1.CameraFocalPoint = C[1]
renderView1.CameraViewUp = C[2]

renderView1.Background = [0.0]*3
ospray = 1
if hasattr(renderView1, 'EnableOSPray'):
    renderView1.EnableOSPRay = ospray
    renderView1.OSPRayRenderer = 'raycaster'
if hasattr(renderView1, 'EnableRayTracing'):
    renderView1.EnableRayTracing = ospray
    renderView1.BackEnd = 'raycaster'
    renderView1.Denoise = 1
renderView1.UseLight = 1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.AmbientSamples = 1
renderView1.SamplesPerPixel = 1 if draft else 10
renderView1.OSPRayMaterialLibrary = materialLibrary1


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

sm = LegacyVTKReader(FileNames=ff)

vf = XDMFReader(FileNames=ffvf)
vf.CellArrayStatus = ['vf']
vf.GridStatus = ['Grid_0']

omm = XDMFReader(FileNames=ffomm)
omm.CellArrayStatus = ['omm']
omm.GridStatus = ['Grid_1']


# list of all sources
vs = [sm, vf, omm]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
sm = ForceTime(sm)
vf = ForceTime(vf)
omm = ForceTime(omm)

# all ForceTime
vft = [sm, vf, omm]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

sm = Calculator(Input=sm)
sm.ResultNormals = 1
sm.AttributeType = 'Point Data'
sm.ResultArrayName = 'normals'
sm.Function = 'nn'

clipbox = Clip(Input=sm)
clipbox.ClipType = 'Scalar'
clipbox.Scalars = ['CELLS', 'cl']
clipbox.Invert = 0
clipbox.Value = -0.1

clipedge = Clip(Input=clipbox)
clipedge.ClipType = 'Box'
clipedge.ClipType.Position = [0.0038, 0.0, 0.0038]
clipedge.ClipType.Length = [1.9924, 1.0, 0.9924]

surface = Clip(Input=clipedge)
surface.ClipType = 'Scalar'
surface.Scalars = ['CELLS', 'cl']
surface.Value = 0.1
surface.Invert = 1
surfaceDisplay = Show(surface, renderView1)
surfaceDisplay.Representation = 'Surface'
surfaceDisplay.ColorArrayName = ['POINTS', '']
surfaceDisplay.Opacity = 1

bubbles = Clip(Input=clipedge)
bubbles.ClipType = 'Scalar'
bubbles.Scalars = ['CELLS', 'cl']
bubbles.Value = 0.1
bubbles.Invert = 0
bubblesDisplay = Show(bubbles, renderView1)
bubblesDisplay.Representation = 'Surface'
bubblesDisplay.ColorArrayName = ['POINTS', '']
bubblesDisplay.AmbientColor = [1.0, 0.7686274509803922, 0.4235294117647059]
bubblesDisplay.DiffuseColor = [1.0, 0.7686274509803922, 0.4235294117647059]


appendAttributes1 = AppendAttributes(Input=[omm, vf])
ommvf = Calculator(Input=appendAttributes1)
ommvf.AttributeType = 'Cell Data'
ommvf.ResultArrayName = 'ommvf'
ommvf.Function = 'omm*(1-vf)^3'
ommvfDisplay = Show(ommvf, renderView1)

ommvfLUT = GetColorTransferFunction('ommvf')
ommvfLUT.AutomaticRescaleRangeMode = 'Never'
ommvfLUT.RGBPoints = [8.0, 0.27058823529411763, 0.4392156862745098, 0.5843137254901961, 70.0, 0.4666666666666667, 0.7294117647058823, 1.0]
ommvfLUT.ColorSpace = 'RGB'
ommvfLUT.NanColor = [1.0, 0.0, 0.0]
ommvfLUT.ScalarRangeInitialized = 1.0
ommvfPWF = GetOpacityTransferFunction('ommvf')
ommvfPWF.Points = [8.0, 0.0, 0.5, 0.0, 70.0, 1.0, 0.5, 0.0]
ommvfPWF.ScalarRangeInitialized = 1

ommvfDisplay.Representation = 'Volume'
ommvfDisplay.AmbientColor = [0.0, 0.0, 0.0]
ommvfDisplay.ColorArrayName = ['CELLS', 'ommvf']
ommvfDisplay.LookupTable = ommvfLUT
ommvfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ommvfDisplay.SelectOrientationVectors = 'None'
ommvfDisplay.OpacityArray = [None, '']
ommvfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ommvfDisplay.DataAxesGrid = 'GridAxesRepresentation'
ommvfDisplay.PolarAxes = 'PolarAxesRepresentation'
ommvfDisplay.ScalarOpacityUnitDistance = 0.03
ommvfDisplay.ScalarOpacityFunction = ommvfPWF
ommvfDisplay.Shade = 1


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
