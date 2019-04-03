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
nfr = len(GetFiles(os.path.join(chdir, "vf_*.xmf")))

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
renderView1.CenterOfRotation = [2.0, 0.5, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [-0.6153985443206067, 2.2621439095968605, 5.021258050312115]
renderView1.CameraFocalPoint = [3.2732823945271643, -0.35788333161888675, -1.7011323242026508]
renderView1.CameraViewUp = [0.18514206388221707, 0.9471305807975574, -0.2620421322984429]
renderView1.CameraParallelScale = 2.1213203435596424
renderView1.Background = [0.32, 0.34, 0.43]
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
def FX(field):
  global chdir
  return F(os.path.join(chdir, "{:}_*.xmf".format(field)))

# Returns reader
def RX(field):
  global gridstatus
# create a new 'CSV Reader'
  ff = FX(field)
  Log("found {:} files with field '{:}'".format(len(ff), field))
  r = XDMFReader(FileNames=ff)
  r.CellArrayStatus = [field]
  r.GridStatus = ['Grid_{:}'.format(gridstatus)]
  gridstatus += 1
  return r

# sources
svf = RX("vf")
svx = RX("vx")

# list of all sources
vs = [svf, svx]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
svf = ForceTime(svf)
svx = ForceTime(svx)

# all ForceTime
vft = [svf, svx]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------

# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(Input=[svf, svx])

# create a new 'Slice'
slice1 = Slice(Input=svx)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [2.0, 0.5, 0.5]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]


# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=svf)

# create a new 'Contour'
contour1 = Contour(Input=cellDatatoPointData1)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1)

# get color transfer function/color map for 'vx'
vxLUT = GetColorTransferFunction('vx')
vxLUT.AutomaticRescaleRangeMode = 'Never'
vxLUT.RGBPoints = [-1.5, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 1.5, 0.705882, 0.0156863, 0.14902]
vxLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'vx']
slice1Display.LookupTable = vxLUT
slice1Display.Opacity = 0.75
slice1Display.Ambient = 0.25


# show data from cellDatatoPointData1
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)
# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Outline'
cellDatatoPointData1Display.ColorArrayName = ['POINTS', '']
cellDatatoPointData1Display.LineWidth = 3.0


# show data from contour1
contour1Display = Show(contour1, renderView1)

contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']

# get color legend/bar for vxLUT in view renderView1
vxLUTColorBar = GetScalarBar(vxLUT, renderView1)
vxLUTColorBar.WindowLocation = 'AnyLocation'
vxLUTColorBar.Position = [0.8683691999132697, 0.0888888888888889]
vxLUTColorBar.Title = 'vx'
vxLUTColorBar.ComponentTitle = ''
vxLUTColorBar.HorizontalTitle = 1
vxLUTColorBar.TitleFontFile = ''
vxLUTColorBar.LabelFontFile = ''
vxLUTColorBar.RangeLabelFormat = '%-#.3g'
vxLUTColorBar.ScalarBarLength = 0.39999999999999997

# set color bar visibility
vxLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'vx'
vxPWF = GetOpacityTransferFunction('vx')
vxPWF.Points = [-1.5, 0.0, 0.5, 0.0, 1.5, 1.0, 0.5, 0.0]
vxPWF.ScalarRangeInitialized = 1


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
