# state file generated using paraview version 5.5.0-RC3

from paraview.simple import *

from glob import glob
import sys
import re

#flog = open("log", 'w')

def Log(s):
  s += "\n"
  #flog.write(s)
  #flog.flush()
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

def GetFiles(pre):
  global base
  p = base + "/{:}_*.xmf".format(pre)
  l = glob(p)
  return natsorted(l)

av = sys.argv
if len(av) < 2:
  sys.stderr.write("usage: {:} basedir\n".format(av[0]))
  exit(1)

# base folder
base = av[1]
# frames to skip
skipfirst = int(av[2]) if len(av) > 2 else 0

# number of frames
nfr = len(GetFiles("vf"))

Log("Using base={:}".format(base))
Log("Skipping first {:} of {:} frames".format(skipfirst, nfr))

if skipfirst >= nfr:
  Log("No frames left")
  exit(0)


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2028, 1186]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [3.14159274101257, 3.14159274101257, 3.14159488677979]
renderView1.StereoType = 0
renderView1.CameraPosition = [20.4106540350961, -4.24286245650427, 9.37676218590444]
renderView1.CameraFocalPoint = [2.31490519108964, 3.59533501902065, 2.08888589162254]
renderView1.CameraViewUp = [-0.317390826200035, 0.139389564802607, 0.937994463026408]
renderView1.CameraParallelScale = 5.44139948298291
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.EnableOSPRay = 1
renderView1.Shadows = 1
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 5
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

def F(pre):
  l = GetFiles(pre)
  assert len(l) == nfr, "found %r files by '%r', expected %r" % (len(l), pre, nfr)
  return l[skipfirst:]

# create a new 'XDMF Reader'
fn = F("vx")
vx_xmf = XDMFReader(FileNames=fn)
vx_xmf.CellArrayStatus = ['vx']
vx_xmf.GridStatus = ['Grid_16805']

# create a new 'XDMF Reader'
fn = F("vy")
vy_xmf = XDMFReader(FileNames=fn)
vy_xmf.CellArrayStatus = ['vy']
vy_xmf.GridStatus = ['Grid_16808']

fn = F("vz")
vz_xmf = XDMFReader(FileNames=fn)
vz_xmf.CellArrayStatus = ['vz']
vz_xmf.GridStatus = ['Grid_16811']

# create a new 'XDMF Reader'
fn = F("vf")
vf_xmf = XDMFReader(FileNames=fn)
vf_xmf.CellArrayStatus = ['vf']
vf_xmf.GridStatus = ['Grid_16802']

Log("Loaded {:} files from {:}".format(len(fn), base))


# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=vf_xmf)

# create a new 'Contour'
contour1 = Contour(Input=cellDatatoPointData1)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from cellDatatoPointData1
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Outline'
#cellDatatoPointData1Display.AmbientColor = [1.0, 1.0, 1.0]
cellDatatoPointData1Display.ColorArrayName = ['POINTS', '']
cellDatatoPointData1Display.LineWidth = 1.0
cellDatatoPointData1Display.OSPRayScaleArray = 'vf'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 0.628318530717958
cellDatatoPointData1Display.SelectScaleArray = 'vf'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'vf'
cellDatatoPointData1Display.GaussianRadius = 0.0314159265358979
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'vf']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'vf']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.SelectionCellLabelFontFile = ''
cellDatatoPointData1Display.SelectionPointLabelFontFile = ''
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
cellDatatoPointData1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.XTitleFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.YTitleFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.ZTitleFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.XLabelFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.YLabelFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
cellDatatoPointData1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.PolarAxes.PolarAxisTitleFontFile = ''
cellDatatoPointData1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.PolarAxes.PolarAxisLabelFontFile = ''
cellDatatoPointData1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.PolarAxes.LastRadialAxisTextFontFile = ''
cellDatatoPointData1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']
contour1Display.OSPRayScaleArray = 'Normals'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.628318548202515
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.GaussianRadius = 0.0314159274101257
contour1Display.SetScaleArray = ['POINTS', 'Normals']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'Normals']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.SelectionCellLabelFontFile = ''
contour1Display.SelectionPointLabelFontFile = ''
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-0.999677538871765, 0.0, 0.5, 0.0, 0.999667048454285, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-0.999677538871765, 0.0, 0.5, 0.0, 0.999667048454285, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
contour1Display.DataAxesGrid.XTitleFontFile = ''
contour1Display.DataAxesGrid.YTitleFontFile = ''
contour1Display.DataAxesGrid.ZTitleFontFile = ''
contour1Display.DataAxesGrid.XLabelFontFile = ''
contour1Display.DataAxesGrid.YLabelFontFile = ''
contour1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour1Display.PolarAxes.PolarAxisTitleFontFile = ''
contour1Display.PolarAxes.PolarAxisLabelFontFile = ''
contour1Display.PolarAxes.LastRadialAxisTextFontFile = ''
contour1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''


# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------

anim = GetAnimationScene()
anim.UpdateAnimationUsingDataTimeSteps()
anim.GoToFirst()

import sys

#nfr = anim.NumberOfFrames
for fr in range(skipfirst,nfr):
  fn = 'aa.{:05d}.png'.format(fr)
  Log("{:}/{:}: {:}".format(fr + 1, nfr, fn))
  SaveScreenshot(fn, renderView1, ImageResolution=[2000,2000])
  anim.GoToNext()

