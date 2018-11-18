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
if len(av) < 3:
  sys.stderr.write('''usage: {:} ch ge
ch: folder with s_*.vtk, interface from ch
ge: folder with u_*.vtk, fields from gerris merged by mfer.cmerge
ba: folder with u_*.vtk, fields from basilisk merged by mfer.cmerge
'''.format(av[0]))
  exit(1)

# base folder
chdir = av[1]  # ch
gedir = av[2]  # ge
badir = av[3]  # ba

# output pattern (:0 substituted by frame number)
bo = "a_{:}.png"

# frames to skip (number of finished frames)
skipfirst = len(glob(bo.format("*")))

# total number of frames
nfr = len(GetFiles(os.path.join(chdir, "s_*.vtk")))
nfr = min(nfr, len(GetFiles(os.path.join(gedir, "u_*.vtk"))))

Log("Using chdir={:} gedir={:}".format(chdir, gedir))
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
renderView1.ViewSize = [1844, 918]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [3.3800271416245082, -0.7474875544264274, 1.6599012405217826]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
renderView1.CameraViewUp = [-0.31739082620003495, 0.13938956480260697, 0.9379944630264078]
renderView1.CameraParallelScale = 0.8660254037844386
renderView1.Background = [0.0, 0.0, 0.0]
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

# create a new 'XDMF Reader'
fn = F(os.path.join(gedir, "u_*.vtk"))
u_00 = LegacyVTKReader(FileNames=fn)

# create a new 'XDMF Reader'
fn = F(os.path.join(badir, "u_*.vtk"))
ub_00 = LegacyVTKReader(FileNames=fn)

# list of all sources
vs = [u_00, s_00, ub_00]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
u_00 = ForceTime(u_00)
s_00 = ForceTime(s_00)
ub_00 = ForceTime(ub_00)

# all ForceTime
vft = [u_00, s_00, ub_00]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------

# create a new 'Transform'
transform1 = Transform(Input=u_00)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [0., 0., 0.]

# create a new 'Contour'
contour1 = Contour(Input=transform1)
contour1.ContourBy = ['POINTS', 'T']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
clpnt = CellDatatoPointData(Input=ub_00)
contour2 = Contour(Input=clpnt)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from transform1
transform1Display = Show(transform1, renderView1)

# trace defaults for the display properties.
transform1Display.Representation = 'Outline'
transform1Display.AmbientColor = [1.0, 1.0, 1.0]
transform1Display.ColorArrayName = [None, '']
transform1Display.LineWidth = 2.0
transform1Display.OSPRayScaleArray = 'P'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 0.1
transform1Display.SelectScaleArray = 'None'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'None'
transform1Display.GaussianRadius = 0.005
transform1Display.SetScaleArray = ['POINTS', 'P']
transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display.OpacityArray = ['POINTS', 'P']
transform1Display.OpacityTransferFunction = 'PiecewiseFunction'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.SelectionCellLabelFontFile = ''
transform1Display.SelectionPointLabelFontFile = ''
transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.ScalarOpacityUnitDistance = 0.054126587736527426

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-1.4124, 0.0, 0.5, 0.0, 28.3541, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-1.4124, 0.0, 0.5, 0.0, 28.3541, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
transform1Display.DataAxesGrid.XTitleFontFile = ''
transform1Display.DataAxesGrid.YTitleFontFile = ''
transform1Display.DataAxesGrid.ZTitleFontFile = ''
transform1Display.DataAxesGrid.XLabelFontFile = ''
transform1Display.DataAxesGrid.YLabelFontFile = ''
transform1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
transform1Display.PolarAxes.PolarAxisTitleFontFile = ''
transform1Display.PolarAxes.PolarAxisLabelFontFile = ''
transform1Display.PolarAxes.LastRadialAxisTextFontFile = ''
transform1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']
contour1Display.DiffuseColor = [1.0, 0.4980392156862745, 0.054901960784313725]
contour1Display.OSPRayScaleArray = 'Normals'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.03803172920884593
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.GaussianRadius = 0.001901586460442296
contour1Display.SetScaleArray = ['POINTS', 'Normals']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'Normals']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.SelectionCellLabelFontFile = ''
contour1Display.SelectionPointLabelFontFile = ''
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-0.9715194702148438, 0.0, 0.5, 0.0, 0.9762166142463684, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-0.9715194702148438, 0.0, 0.5, 0.0, 0.9762166142463684, 1.0, 0.5, 0.0]

# show data from contour2
contour2Display = Show(contour2, renderView1)

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = [None, '']
contour2Display.DiffuseColor = [0., 0.8, 0.]
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


# show data from s_00
s_00Display = Show(s_00, renderView1)

# trace defaults for the display properties.
s_00Display.Representation = 'Surface'
s_00Display.ColorArrayName = ['POINTS', '']
s_00Display.DiffuseColor = [0.12156862745098039, 0.4666666666666667, 0.7058823529411765]
s_00Display.Opacity = 0.75
s_00Display.Ambient = 0.4
s_00Display.OSPRayScaleFunction = 'PiecewiseFunction'
s_00Display.SelectOrientationVectors = 'None'
s_00Display.ScaleFactor = 0.015503498911857606
s_00Display.SelectScaleArray = 'c'
s_00Display.GlyphType = 'Arrow'
s_00Display.GlyphTableIndexArray = 'c'
s_00Display.GaussianRadius = 0.0007751749455928803
s_00Display.SetScaleArray = [None, '']
s_00Display.ScaleTransferFunction = 'PiecewiseFunction'
s_00Display.OpacityArray = [None, '']
s_00Display.OpacityTransferFunction = 'PiecewiseFunction'
s_00Display.DataAxesGrid = 'GridAxesRepresentation'
s_00Display.SelectionCellLabelFontFile = ''
s_00Display.SelectionPointLabelFontFile = ''
s_00Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_00Display.DataAxesGrid.XTitleFontFile = ''
s_00Display.DataAxesGrid.YTitleFontFile = ''
s_00Display.DataAxesGrid.ZTitleFontFile = ''
s_00Display.DataAxesGrid.XLabelFontFile = ''
s_00Display.DataAxesGrid.YLabelFontFile = ''
s_00Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_00Display.PolarAxes.PolarAxisTitleFontFile = ''
s_00Display.PolarAxes.PolarAxisLabelFontFile = ''
s_00Display.PolarAxes.LastRadialAxisTextFontFile = ''
s_00Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

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
