#!/usr/bin/env pvbatch

# Plots interfaces from ch and ge.
# $1: folder with s_*.vtk, interface from ch
# $2: folder with u_*.vtk, fields from gerris, merged by mfer.cmerge
# Requires pvbatch from Paraview 5.5.1.
# Output:
# a_*.png in current folder

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


# Sets time of datasets to step i
def SetTime(i):
    global vft, vt
    for j in range(len(vft)):
        s = vft[j]
        s.ForcedTime = vt[j][i]
        s.UpdatePipeline()


av = sys.argv
if len(av) < 2 or av[1] == '-h':
    sys.stderr.write('''usage: {:} [s_*.vtk]
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

# s input
ff = natsorted(av[1:])
# s basename
ffb = list(map(os.path.basename, ff))
# s dirname
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
# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000, 1000]
renderView1.OrientationAxesVisibility = 0

renderView1.CameraPosition = [1.7536929606099154, 3.44405411327432, 1.4829402975071122]
renderView1.CameraFocalPoint = [0.47975814772145614, 0.5400876539226072, 0.4151288356011742]
renderView1.CameraViewUp = [-0.11344083319471632, -0.2986712073425839, 0.9475899362427884]
renderView1.CameraParallelScale = 0.8660254037844386

#renderView1.CameraPosition = [3.38060827699219, 2.283201516563153, 1.925719089019077]
#renderView1.CameraFocalPoint = [0.5000000000000002, 0.5000000000000016, 0.9987755119800572]
#renderView1.CameraViewUp = [-0.22201257264706833, -0.14274411445459043, 0.9645385090162057]
#renderView1.CameraParallelScale = 0.7513019528597746

renderView1.Background = [1.0, 1.0, 1.0]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# list of all sources
vs = []
# list of timesteps arrays
vt = []
# list of ForceTime sources
vft = []

def A(s):
    global vs, vt, vft
    if s is not None:
        vs.append(s)
        vt.append(np.array(s.TimestepValues))
        s = ForceTime(s)
        vft.append(s)
    return s

# create a new 'CSV Reader'
s_00 = LegacyVTKReader(FileNames=ff)
s_00 = A(s_00)

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------

# show data from s_00
s_00Display = Show(s_00, renderView1)

# trace defaults for the display properties.
s_00Display.Representation = 'Surface'
#s_00Display.Representation = 'Surface'
s_00Display.AmbientColor = [
        0.12156862745098039, 0.4666666666666667, 0.7058823529411765]
s_00Display.ColorArrayName = ['POINTS', '']
s_00Display.LineWidth = 3.0

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------

#####################################################
### END OF STATE FILE
#####################################################

anim = GetAnimationScene()
anim.UpdateAnimationUsingDataTimeSteps()
anim.GoToFirst()

# XXX: workaround for blank contour on first frame
anim.GoToNext()
anim.GoToPrevious()

for i in list(range(len(ss))):
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue

    SetTime(i)

    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1)

exit(0)
