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

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000,1000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 1.5, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition =\
[1.326166076499316, 1.581412512627694, 3.5052432288262967]
renderView1.CameraFocalPoint =\
[0.45722849374168073, 0.4442839115593683, 0.4806707888589165]
renderView1.CameraViewUp =\
[-0.08433283755119116, 0.9404371241304461, -0.3293417496580418]
renderView1.CameraParallelScale = 1.6583123951777




renderView1.Background = [0.0, 0.0, 0.0]
renderView1.EnableOSPRay = ospray
renderView1.Shadows = 1
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 10
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

surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = 'Surface'
surfDisplay.AmbientColor = [0.0, 0.0, 0.0]
surfDisplay.Opacity = 0.5
surfDisplay.ColorArrayName = [None, '']

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
