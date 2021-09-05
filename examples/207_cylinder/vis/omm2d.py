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
# a_*.png in current directory
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

par = {}
exec(open("../par.py").read(), None, par)
extent = par["extent"]
lx = extent
ly = extent * 0.5
radius = 0.5
center = [lx * 0.3, ly * 0.5, 0]

C1 = [
    [lx * 0.3, ly * 0.5, 10],
    [lx * 0.3, ly * 0.5, 0],
    [0, 1, 0],
    2,
    ]

CC = [C1]
C = CC[cam - 1]

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2000, 1000]
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
renderView1.UseLight = 0
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.OSPRayMaterialLibrary = materialLibrary1


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

omm = XDMFReader(FileNames=ff)
ommname = 'omz'
omm.CellArrayStatus = [ommname]
omm.GridStatus = ['Grid_1']


# list of all sources
vs = [omm]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
omm = ForceTime(omm)

# all ForceTime
vft = [omm]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

k = 0.5
ommDisplay = Show(omm, renderView1)

ommLUT = GetColorTransferFunction(ommname)
ommLUT.AutomaticRescaleRangeMode = 'Never'
if ommname == 'omz':
    ommLUT.RGBPoints = [-25.0*k, 0.231373, 0.298039, 0.752941, 0*k, 0.865003, 0.865003, 0.865003, 25.0*k, 0.705882, 0.0156863, 0.14902]
else:
    ommLUT.RGBPoints = [0.0*k, 0.231373, 0.298039, 0.752941, 12.500000000000004*k, 0.865003, 0.865003, 0.865003, 25.0*k, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0

ommPWF = GetOpacityTransferFunction(ommname)
if ommname == 'omz':
    ommPWF.Points = [-25.0*k, 1.0, 0.5, 0.0, 25.0*k, 1.0, 0.5, 0.0]
else:
    ommPWF.Points = [0.0*k, 1.0, 0.5, 0.0, 25.0*k, 1.0, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
ommDisplay.Representation = 'Surface'
ommDisplay.ColorArrayName = ['CELLS', ommname]
ommDisplay.LookupTable = ommLUT

ebvtk = LegacyVTKReader(FileNames=['../eb.vtk'])
ebvtk = Clip(Input=ebvtk)
ebvtk.ClipType = 'Scalar'
ebvtk.Scalars = ['CELLS', 'face']
ebvtk.Value = 0.5
ebvtkDisplay = Show(ebvtk, renderView1)
ebvtkDisplay.Representation = 'Wireframe'
ebvtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
ebvtkDisplay.ColorArrayName = ['POINTS', '']
ebvtkDisplay.DiffuseColor = [0.0, 0.0, 0.0]
ebvtkDisplay.LineWidth = 5.

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
