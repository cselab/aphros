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
    sys.stderr.write('''usage: {:} [vf_*.xmf]
Plots bubbles with vorticity.
Current folder:
omm_*.xmf
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
# omm input
ffomm = ["omm_{:04d}.xmf".format(s) for s in ss]

# append dirname
for i in range(len(ss)):
    ffomm[i] = os.path.join(ffd[i], ffomm[i])

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

renderView1.CenterOfRotation = \
    [2.0002396481577307, 0.49982515879673883, 0.5001013136352412]
renderView1.CameraPosition = \
    [4.1725221459017545, 4.290558691659238, 2.1117181099306443]
renderView1.CameraFocalPoint = \
    [0.37576486475267523, -2.37787637240639, -0.7371351849693569]
renderView1.CameraViewUp = \
    [-0.17297883552233423, -0.3020165314331066, 0.9374776462414731]
renderView1.CameraParallelScale = 1.1625484322433086
renderView1.CameraParallelProjection = 1

renderView1.Shadows = 0
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableOSPRay = 1
renderView1.AmbientSamples = 10
renderView1.SamplesPerPixel = 10
renderView1.LightScale = 0.8
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

plane1 = Plane()
plane1.Origin = [0.0, 0.0, 0.0]
plane1.Point1 = [0.0, 1.0, 0.0]
plane1.Point2 = [4.0, 0.0, 0.0]
plane1Display = Show(plane1, renderView1)
plane1Display.Representation = 'Surface'
plane1Display.ColorArrayName = [None, '']

plane2 = Plane()
plane2.Origin = [0.0, 0.0, 1.0]
plane2.Point1 = [0.0, 1.0, 1.0]
plane2.Point2 = [4.0, 0.0, 1.0]
plane2Display = Show(plane2, renderView1)
plane2Display.Representation = 'Surface'
plane2Display.ColorArrayName = [None, '']
plane2Display.Opacity = 0.5

plane3 = Plane()
plane3.Origin = [0.0, 0.0, 0.0]
plane3.Point1 = [4.0, 0.0, 0.0]
plane3.Point2 = [0.0, 0.0, 1.0]
plane3Display = Show(plane3, renderView1)
plane3Display.Representation = 'Surface'
plane3Display.ColorArrayName = [None, '']

surfshow = Show(surf, renderView1)
surfshow.Representation = 'Surface'
surfshow.ColorArrayName = ['CELLS', '']


#####################################################
### END OF STATE FILE
#####################################################

SetTime(1)
SaveScreenshot("tmp.png", renderView1)

for i in list(range(len(ss))):
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue

    SetTime(i)

    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1)

exit(0)
