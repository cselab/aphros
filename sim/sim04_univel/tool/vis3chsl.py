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
renderView1.ViewSize = [1000, 2000]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.StereoType = 0
renderView1.CameraPosition = [0.25, -10, 0.5]
renderView1.CameraFocalPoint = [0.25, 0.5, 0.5]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.5
renderView1.CameraParallelProjection = 1
renderView1.Background = [1, 1, 1]
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
svy = RX("vy")
svz = RX("vz")
sp = RX("p")

# list of all sources
vs = [svf, svx, svy, svz, sp]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
svf = ForceTime(svf)
svx = ForceTime(svx)
svy = ForceTime(svy)
svz = ForceTime(svz)
sp = ForceTime(sp)

# all ForceTime
vft = [svf, svx, svy, svz, sp]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------

# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(Input=[sp, svf, svx, svy, svz])

# create a new 'Calculator'
calculator1 = Calculator(Input=appendAttributes1)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'vel'
# XXX: vy=0
calculator1.Function = '(iHat*vx+jHat*vy*0+kHat*vz)'

calculator2 = Calculator(Input=calculator1)
calculator2.AttributeType = 'Cell Data'
calculator2.ResultArrayName = 'vm'
calculator2.Function = 'mag(vel)'

# create a new 'Calculator'
calculator3 = Calculator(Input=calculator2)
calculator3.AttributeType = 'Cell Data'
calculator3.ResultArrayName = 'vmm'
calculator3.Function = 'vm*(1-vf)'

# create a new 'Slice'
slice1 = Slice(Input=calculator3)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.25, 0.5, 0.5]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=slice1)

# create a new 'Contour'
contour1 = Contour(Input=cellDatatoPointData1)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Glyph'
glyph1 = Glyph(Input=slice1, GlyphType='Arrow')
glyph1.OrientationArray = ['CELLS', 'vel']
glyph1.ScaleArray = ['CELLS', 'vmm']
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'Every Nth Point'
glyph1.ScaleFactor = 0.4
glyph1.Stride = 99

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.AutomaticRescaleRangeMode = 'Never'
pLUT.RGBPoints = [-0.05, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 0.05, 0.705882, 0.0156863, 0.14902]
pLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.AmbientColor = [0.0, 0.0, 0.0]
slice1Display.ColorArrayName = ['CELLS', 'p']
slice1Display.LookupTable = pLUT

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.ColorArrayName = [None, '']
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]
contour1Display.PointSize = 30.0
contour1Display.LineWidth = 10.0
contour1Display.RenderPointsAsSpheres = 1
contour1Display.OSPRayScaleArray = 'p'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'vel'
contour1Display.ScaleFactor = 0.0089234858751297
contour1Display.GlyphType = 'Arrow'

# show data from glyph1
glyph1Display = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
glyph1Display.ColorArrayName = ['POINTS', '']
glyph1Display.DiffuseColor = [0.0, 0.0, 0.0]


# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')
pPWF.Points = [-0.05, 0.0, 0.5, 0.0, 0.05, 1.0, 0.5, 0.0]
pPWF.ScalarRangeInitialized = 1

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
