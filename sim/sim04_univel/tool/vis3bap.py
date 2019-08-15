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
  sys.stderr.write('''usage: {:} ba
ba: folder with u_*.vtk, fields from basilisk merged by mfer.cmerge
'''.format(av[0]))
  exit(1)

# base folder
badir = av[1]  # ba

# output pattern (:0 substituted by frame number)
bo = "a_{:}.png"

# frames to skip (number of finished frames)
skipfirst = len(glob(bo.format("*")))

# total number of frames
nfr = len(GetFiles(os.path.join(badir, "u_*.vtk")))

Log("Using badir={:}".format(badir))
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
renderView1.ViewSize = [800, 800]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.5, -3.8460652149512318, 0.5]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
renderView1.CameraViewUp = [0., 0., 1.]
renderView1.CameraParallelScale = 0.499939956838161
renderView1.CameraParallelProjection = 1
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
# create a new 'XDMF Reader'
fn = F(os.path.join(badir, "u_*.vtk"))
ub_00 = LegacyVTKReader(FileNames=fn)

# list of all sources
vs = [ub_00]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
ub_00 = ForceTime(ub_00)

# all ForceTime
vft = [ub_00]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------


# create a new 'Contour'
clpnt = CellDatatoPointData(Input=ub_00)
contour2 = Contour(Input=clpnt)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------


# show data from contour2
contour2Display = Show(contour2, renderView1)

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = [None, '']
contour2Display.DiffuseColor = [0.12156862745098039, 0.4666666666666667, 0.7058823529411765]
contour2Display.OSPRayScaleArray = 'Normals'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.03803172920884593
contour2Display.SelectScaleArray = 'None'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'None'
contour2Display.GaussianRadius = 0.001901586460442296
contour2Display.SetScaleArray = ['POINTS', 'Normals']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'Normals']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.SelectionCellLabelFontFile = ''
contour2Display.SelectionPointLabelFontFile = ''
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [-0.9715194702148438, 0.0, 0.5, 0.0, 0.9762166142463684, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [-0.9715194702148438, 0.0, 0.5, 0.0, 0.9762166142463684, 1.0, 0.5, 0.0]



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
  SaveScreenshot(fn, renderView1, ImageResolution=[1000,1000])
  anim.GoToNext()

exit(0)
