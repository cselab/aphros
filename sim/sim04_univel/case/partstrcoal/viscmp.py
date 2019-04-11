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
  sys.stderr.write('''usage: {:} ch1 ch2
Plots difference ch2-ch1
ch1,ch2: folders with s_*.vtk, interface from ch
'''.format(av[0]))
  exit(1)

# base folder
chdir = av[1]  # ch
chdir2 = av[2]  # ch

# output pattern (:0 substituted by frame number)
bo = "a_{:}.png"

# frames to skip (number of finished frames)
skipfirst = len(glob(bo.format("*")))

# total number of frames
nfr = len(GetFiles(os.path.join(chdir, "vf_*.xmf")))
nfr2 = len(GetFiles(os.path.join(chdir, "vf_*.xmf")))
nfr = min(nfr, nfr2)

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
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.25, 0.25, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [-0.1894920791036435, -0.8468369066649195, 0.28994809267666166]
renderView1.CameraFocalPoint = [0.6774136004022462, 1.3394285117826572, 0.03150220996103598]
renderView1.CameraViewUp = [0.04055106344301122, 0.10142659343302993, 0.9940162259230036]
renderView1.CameraParallelScale = 0.6123724356957945
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

gridstatus = 0

# Returns list of files by glob pattern p skipping skipfirst
def F(p):
  global skipfirst
  l = GetFiles(p)
  #assert len(l) >= nfr, "found %r files by '%r', expected at least %r" % (len(l), p, nfr)
  return l[skipfirst:]

# Returns list of files by glob pattern p skipping skipfirst
def FX(field, d):
  return F(os.path.join(d, "{:}_*.xmf".format(field)))

# Returns list of files by glob pattern p skipping skipfirst
def FV(field, d):
  return F(os.path.join(d, "{:}_*.vtk".format(field)))

# Returns reader
def RX(field, d):
  global gridstatus
# create a new 'CSV Reader'
  ff = FX(field, d)
  Log("found {:} files with field '{:}'".format(len(ff), field))
  r = XDMFReader(FileNames=ff)
  r.CellArrayStatus = [field]
  r.GridStatus = ['Grid_{:}'.format(gridstatus)]
  gridstatus += 1
  return r

# Returns reader
def RV(field, d):
  ff = FV(field, d)
  Log("found {:} files with field '{:}'".format(len(ff), field))
  r = LegacyVTKReader(FileNames=ff)
  return r

# sources
field = "vz"
sv1 = RX(field, chdir)
sv2 = RX(field, chdir2)
ss1 = RV("s", chdir)
ss2 = RV("s", chdir2)

# list of all sources
vs = [sv1, sv2, ss1, ss2]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
sv1 = ForceTime(sv1)
sv2 = ForceTime(sv2)
ss1 = ForceTime(ss1)
ss2 = ForceTime(ss2)

# all ForceTime
vft = [sv1, sv2, ss1, ss2]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------

calcv1 = Calculator(Input=sv1)
calcv1.AttributeType = 'Cell Data'
calcv1.ResultArrayName = 'v1'
calcv1.Function = field

calcv2 = Calculator(Input=sv2)
calcv2.AttributeType = 'Cell Data'
calcv2.ResultArrayName = 'v2'
calcv2.Function = field

# create a new 'Append Attributes'
appnd = AppendAttributes(Input=[calcv1, calcv2])

calcdv = Calculator(Input=appnd)
calcdv.AttributeType = 'Cell Data'
calcdv.ResultArrayName = 'dv'
calcdv.Function = 'v2-v1'

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=calcdv)
resampleToImage1.SamplingDimensions = [129, 129, 257]
resampleToImage1.SamplingBounds = [0.0, 0.5, 0.0, 0.5, 0.0, 1.0]


# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
resampleToImage1Display = Show(resampleToImage1, renderView1)

# get color transfer function/color map for 'dvz'
dvzLUT = GetColorTransferFunction('dv')
dvzLUT.AutomaticRescaleRangeMode = 'Never'
dvzLUT.RGBPoints = [-0.05, 0.231373, 0.298039, 0.752941, 1.3877787807814457e-17, 0.865003, 0.865003, 0.865003, 0.05, 0.705882, 0.0156863, 0.14902]
dvzLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'dvz'
dvzPWF = GetOpacityTransferFunction('dv')
dvzPWF.Points = [-0.05, 1.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.05, 1.0, 0.5, 0.0]
dvzPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
resampleToImage1Display.Representation = 'Volume'
resampleToImage1Display.AmbientColor = [0.0, 0.0, 0.0]
resampleToImage1Display.ColorArrayName = ['POINTS', 'dv']
resampleToImage1Display.LookupTable = dvzLUT
resampleToImage1Display.OSPRayScaleArray = 'dv'
resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
resampleToImage1Display.SelectOrientationVectors = 'None'
resampleToImage1Display.ScaleFactor = 0.1
resampleToImage1Display.SelectScaleArray = 'None'
resampleToImage1Display.GlyphType = 'Arrow'
resampleToImage1Display.GlyphTableIndexArray = 'None'
resampleToImage1Display.GaussianRadius = 0.005
resampleToImage1Display.SetScaleArray = ['POINTS', 'dv']
resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.OpacityArray = ['POINTS', 'dv']
resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleToImage1Display.ScalarOpacityFunction = dvzPWF
# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
resampleToImage1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
resampleToImage1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
resampleToImage1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from s_00
s_00Display = Show(ss1, renderView1)

# trace defaults for the display properties.
s_00Display.Representation = 'Surface'
s_00Display.AmbientColor = [0.0, 0.0, 0.0]
s_00Display.ColorArrayName = ['POINTS', '']
s_00Display.Opacity = 0.5

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
