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
    sys.stderr.write('''usage: {:} [sm_*.vtk]
Plots isosurface.
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

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.OrientationAxesVisibility = 0
renderView1.CameraPosition =\
[2.18498889205859, 1.6310649466385199, 2.5118462518730618]
renderView1.CameraFocalPoint =\
[0.11174881869864335, -0.5360846098832542, -1.1484055486585474]
renderView1.CameraViewUp =\
[-0.27081993085019135, 0.8877545682912452, -0.3722208907790139]
renderView1.CameraParallelScale = 1.224744871391589
renderView1.CameraParallelProjection = 0
renderView1.Background = [1.0]*3
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
surf.ResultArrayName = 'normals'
surf.Function = 'nn'

water = Clip(Input=surf)
water.ClipType = 'Scalar'
water.Scalars = ['CELLS', 'cl']
water.Value = 1.0
water.Invert = 1
waterDisplay = Show(water, renderView1)
waterDisplay.Representation = 'Surface'
waterDisplay.ColorArrayName = [None, '']
waterDisplay.Opacity = 0.5

bubbles = Clip(Input=surf)
bubbles.ClipType = 'Scalar'
bubbles.Scalars = ['CELLS', 'cl']
bubbles.Value = 1.0
bubbles.Invert = 0
bubblesDisplay = Show(bubbles, renderView1)
bubblesDisplay.Representation = 'Surface'
bubblesDisplay.ColorArrayName = [None, '']
bubblesDisplay.DiffuseColor = [1.0, 0.0, 0.0]

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
