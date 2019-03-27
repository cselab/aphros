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
# Returns sorted list of files in base by pattern pre_*.xmf

# Returns sorted list of files by glob pattern p
def GetFiles(p):
  l = glob(p)
  return natsorted(l)

av = sys.argv
if len(av) < 2:
  sys.stderr.write('''usage: {:} ch
ch: folder with s_*.vtk, interface from ch
'''.format(av[0]))
  exit(1)

# base folder
chdir = av[1]  # ch

# output pattern (:0 substituted by frame number)
bo = "a_{:}.png"

# frames to skip (number of finished frames)
skipfirst = len(glob(bo.format("*")))

# total number of frames
nfr = len(GetFiles(os.path.join(chdir, "s_*.vtk")))

Log("Using chdir={:}".format(chdir))
Log("Skipping first {:} of {:} frames".format(skipfirst, nfr))

if skipfirst >= nfr:
  Log("No frames left")
  exit(0)

#####################################################
### BEGIN OF STATE FILE
#####################################################

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1280, 720]
renderView1.CameraPosition = [-0.6815454304186566, 1.1535098106881623, 3.2532374959707173]
renderView1.CameraFocalPoint = [0.01636127691530816, 0.812988767813245, 2.066255492031417]
renderView1.CameraViewUp = [0.1354861087701643, 0.9706311712039707, -0.19879296722351014]
renderView1.CameraViewAngle = 26.885119691299796
renderView1.CameraParallelScale = 0.8445041384116267
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.OSPRayMaterialLibrary = materialLibrary1


# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# Returns list of files by glob pattern p skipping skipfirst
def F(p):
  global skipfirst
  l = GetFiles(p)
  #assert len(l) >= nfr, "found %r files by '%r', expected at least %r" % (len(l), p, nfr)
  return l[skipfirst:]

# create a new 'CSV Reader'
fn = F(os.path.join(chdir, "s_*.vtk"))
s_00 = LegacyVTKReader(FileNames=fn)

# list of all sources
vs = [s_00]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
s_00 = ForceTime(s_00)

# all ForceTime
vft = [s_00]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------


# show data from s_00
s_00Display = Show(s_00, renderView1)

# trace defaults for the display properties.
s_00Display.Representation = 'Surface'
s_00Display.ColorArrayName = ['POINTS', '']
s_00Display.DiffuseColor = [0.1, 0.5, 0.7]
s_00Display.Opacity = 1.
s_00Display.Ambient = 0.

box1 = Box()
box1.XLength = 2.0
box1.Center = [1., 0.5, 0.5]
box1Display = Show(box1, renderView1)
box1Display.Representation = 'Wireframe'
box1Display.AmbientColor = [0.0, 0.0, 0.0]
box1Display.LineWidth = 2.0


#####################################################
### END OF STATE FILE
#####################################################

anim = GetAnimationScene()
anim.UpdateAnimationUsingDataTimeSteps()
anim.GoToFirst()

# XXX: workaround for blank contour on first frame
anim.GoToNext()
anim.GoToPrevious()

import sys

for fr in range(skipfirst,nfr):
  dfr = fr - skipfirst
  for k in range(len(vft)):
    vft[k].SetPropertyWithName("ForcedTime", vt[k][dfr])
  fn = bo.format("{:04d}".format(fr))
  Log("{:}/{:}: {:}".format(fr + 1, nfr, fn))
  SaveScreenshot(fn, renderView1)
  anim.GoToNext()

exit(0)
