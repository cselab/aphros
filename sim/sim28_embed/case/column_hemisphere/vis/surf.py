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

ext = 1.35
# view
C1 = [
[0.5000021234154701*ext, 1*ext, 4.26489794514011*ext],
[0.5000021234154701*ext, 1*ext, -4.6437767967976145*ext],
[0, 1, 0],
1,
    ]


CC = [C1]
C = CC[cam - 1]

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [500, 1000]
renderView1.OrientationAxesVisibility = 0
renderView1.CameraPosition = C[0]
renderView1.CameraFocalPoint = C[1]
renderView1.CameraViewUp = C[2]
renderView1.Background = [1.0]*3
renderView1.OSPRayMaterialLibrary = materialLibrary1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.CameraParallelScale = C[3]
renderView1.CameraParallelProjection = 0


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

def rgb(r, g, b):
    m = 255.
    return [r/m, g/m, b/m]

bubblesDisplay = Show(surf, renderView1)
bubblesDisplay.Representation = 'Surface'
bubblesDisplay.ColorArrayName = [None, '']
bubblesDisplay.DiffuseColor = rgb(255, 127, 14)
bubblesDisplay.Ambient = 0.2

ebvtk = LegacyVTKReader(FileNames=['../eb.vtk'])
ebvtk = Clip(Input=ebvtk)
ebvtk.ClipType = 'Scalar'
ebvtk.Scalars = ['CELLS', 'face']
ebvtk.Value = 0.5

ebvtk = Clip(Input=ebvtk)
ebvtk.ClipType = 'Plane'
ebvtk.ClipType.Origin = [0.5, 0.5, 0.75]
ebvtk.ClipType.Normal = [0.0, 0.0, 1.0]

ebvtkDisplay = Show(ebvtk, renderView1)
ebvtkDisplay.Representation = 'Surface'
ebvtkDisplay.AmbientColor = [1.0, 1.0, 1.0]
ebvtkDisplay.ColorArrayName = ['POINTS', '']
ebvtkDisplay.DiffuseColor = [1.0, 1.0, 1.0]

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
