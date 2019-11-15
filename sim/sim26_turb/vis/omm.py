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

draft = CheckFlag('-draft')

# sm input
ff = natsorted(av[1:])
# sm basename
ffb = list(map(os.path.basename, ff))
# sm dirname
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

materialLibrary1 = GetMaterialLibrary()


# standard view from corner
C1 = [
     [6.43080086005585, -15.166333096764594, 9.914876482661194],
     [3.0001989267510476, 4.2895772779323025, 2.724262644294233],
     [-0.05939117461388465, 0.33682408883346476, 0.9396926207859086],
     4.2537082476917885,
    ]

CC = [C1]
C = CC[cam - 1]

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1800, 2000]
if draft:
    renderView1.ViewSize[0] = renderView1.ViewSize[0] // 2
    renderView1.ViewSize[1] = renderView1.ViewSize[1] // 2
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 0
renderView1.CameraPosition = C[0]
renderView1.CameraFocalPoint = C[1]
renderView1.CameraViewUp = C[2]
renderView1.CameraParallelScale = C[3]
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.]*3
renderView1.UseLight = 1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.OSPRayMaterialLibrary = materialLibrary1


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

sm = LegacyVTKReader(FileNames=ff)

omm = XDMFReader(FileNames=ffomm)
omm.CellArrayStatus = ['omm']
omm.GridStatus = ['Grid_1']


# list of all sources
vs = [sm, omm]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
sm = ForceTime(sm)
omm = ForceTime(omm)

# all ForceTime
vft = [sm, omm]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

sm = Calculator(Input=sm)
sm.ResultNormals = 1
sm.AttributeType = 'Point Data'
sm.ResultArrayName = 'normals'
sm.Function = 'nn'

surface = sm
surfaceDisplay = Show(surface, renderView1)
surfaceDisplay.Representation = 'Surface'
surfaceDisplay.ColorArrayName = ['CELLS', '']
surfaceDisplay.Opacity = 1

k=10
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [1.*k, 0.231373, 0.298039, 0.752941, 3.*k, 0.865003, 0.865003, 0.865003, 5.*k, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [1.0*k, 0.0, 0.5, 0.0, 3.*k, 0.07236842066049576, 0.5, 0.0, 5.*k, 0.2565789520740509, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1

gradDisplay = Show(omm, renderView1)
gradDisplay.Representation = 'Volume'
gradDisplay.ColorArrayName = ['CELLS', 'omm']
gradDisplay.LookupTable = ommLUT
gradDisplay.ScalarOpacityUnitDistance = 0.04251092259923938
gradDisplay.ScalarOpacityFunction = ommPWF

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
