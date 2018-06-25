# state file generated using paraview version 5.5.0-RC3

from paraview.simple import *

from glob import glob
import sys
import re

flog = open("log", 'w')

def Log(s):
  s += "\n"
  flog.write(s)
  flog.flush()
  o = sys.stdout
  o.write(s)
  o.flush()

# natural sort
def natkey(s, _nsre=re.compile('([0-9]+)')):
      return [int(text) if text.isdigit() else text.lower()
                      for text in re.split(_nsre, s)]

def natsorted(v):
  return sorted(v, key=natkey)


# read data folder
av = sys.argv
if len(av) < 2:
  sys.stderr.write("usage: {:} basedir\n".format(av[0]))
  exit(1)
else:
  base = av[1]

skipfirst = 0
if len(av) > 2:
  skipfirst = int(av[2])

Log("Using base={:}".format(base))
Log("skipping first {:} files".format(skipfirst))


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

# create a new 'XDMF Reader'

# Returns sorted list of files found by pattern pre_*.xmf
# pre: prefix
def F(pre):
  p = base + "/{:}_*.xmf".format(pre)
  l = glob(p)
  assert len(l), "no files found by pattern %r" % p
  Log("Found {:} files by pattern '{:}'".format(len(l), p))
  return natsorted(l)[skipfirst:]

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


# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(Input=[vx_xmf, vy_xmf, vz_xmf])

# create a new 'Cell Data to Point Data'
cellDatatoPointData2 = CellDatatoPointData(Input=appendAttributes1)

# create a new 'Calculator'
calculator1 = Calculator(Input=cellDatatoPointData2)
calculator1.ResultArrayName = 'vel'
calculator1.Function = 'vx*iHat+vy*jHat+vz*kHat'

# create a new 'Gradient Of Unstructured DataSet'
gradientOfUnstructuredDataSet1 = GradientOfUnstructuredDataSet(Input=calculator1)
gradientOfUnstructuredDataSet1.ScalarArray = ['POINTS', 'vel']
gradientOfUnstructuredDataSet1.ComputeGradient = 0
gradientOfUnstructuredDataSet1.ComputeVorticity = 1

# create a new 'Extract Component'
extractComponent1 = ExtractComponent(Input=gradientOfUnstructuredDataSet1)
extractComponent1.InputArray = ['POINTS', 'Vorticity']
extractComponent1.Component = 2
extractComponent1.OutputArrayName = 'omz'

# create a new 'Contour'
contour2 = Contour(Input=extractComponent1)
contour2.ContourBy = ['POINTS', 'omz']
contour2.Isosurfaces = [2.0]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour3 = Contour(Input=extractComponent1)
contour3.ContourBy = ['POINTS', 'omz']
contour3.Isosurfaces = [-2.0]
contour3.PointMergeMethod = 'Uniform Binning'

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

# show data from contour2
contour2Display = Show(contour2, renderView1)

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', '']
contour2Display.DiffuseColor = [1.0, 0.25882352941176473, 0.3568627450980392]
contour2Display.Opacity = 0.5
contour2Display.OSPRayScaleArray = 'vx'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'vel'
contour2Display.ScaleFactor = 0.628318548202515
contour2Display.SelectScaleArray = 'vx'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'vx'
contour2Display.GaussianRadius = 0.0314159274101257
contour2Display.SetScaleArray = ['POINTS', 'vx']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'vx']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.SelectionCellLabelFontFile = ''
contour2Display.SelectionPointLabelFontFile = ''
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [-1.5043346601308554, 0.0, 0.5, 0.0, 1.9631953650664902, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [-1.5043346601308554, 0.0, 0.5, 0.0, 1.9631953650664902, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
contour2Display.DataAxesGrid.XTitleFontFile = ''
contour2Display.DataAxesGrid.YTitleFontFile = ''
contour2Display.DataAxesGrid.ZTitleFontFile = ''
contour2Display.DataAxesGrid.XLabelFontFile = ''
contour2Display.DataAxesGrid.YLabelFontFile = ''
contour2Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour2Display.PolarAxes.PolarAxisTitleFontFile = ''
contour2Display.PolarAxes.PolarAxisLabelFontFile = ''
contour2Display.PolarAxes.LastRadialAxisTextFontFile = ''
contour2Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from contour3
contour3Display = Show(contour3, renderView1)

# trace defaults for the display properties.
contour3Display.Representation = 'Surface'
contour3Display.ColorArrayName = ['POINTS', '']
contour3Display.DiffuseColor = [0.48627450980392156, 0.5568627450980392, 1.0]
contour3Display.Opacity = 0.5
contour3Display.OSPRayScaleArray = 'vx'
contour3Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour3Display.SelectOrientationVectors = 'vel'
contour3Display.ScaleFactor = 0.6283185482025146
contour3Display.SelectScaleArray = 'vx'
contour3Display.GlyphType = 'Arrow'
contour3Display.GlyphTableIndexArray = 'vx'
contour3Display.GaussianRadius = 0.031415927410125735
contour3Display.SetScaleArray = ['POINTS', 'vx']
contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
contour3Display.OpacityArray = ['POINTS', 'vx']
contour3Display.OpacityTransferFunction = 'PiecewiseFunction'
contour3Display.DataAxesGrid = 'GridAxesRepresentation'
contour3Display.SelectionCellLabelFontFile = ''
contour3Display.SelectionPointLabelFontFile = ''
contour3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour3Display.ScaleTransferFunction.Points = [-1.3271294794117883, 0.0, 0.5, 0.0, 1.5760667770943608, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour3Display.OpacityTransferFunction.Points = [-1.3271294794117883, 0.0, 0.5, 0.0, 1.5760667770943608, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
contour3Display.DataAxesGrid.XTitleFontFile = ''
contour3Display.DataAxesGrid.YTitleFontFile = ''
contour3Display.DataAxesGrid.ZTitleFontFile = ''
contour3Display.DataAxesGrid.XLabelFontFile = ''
contour3Display.DataAxesGrid.YLabelFontFile = ''
contour3Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour3Display.PolarAxes.PolarAxisTitleFontFile = ''
contour3Display.PolarAxes.PolarAxisLabelFontFile = ''
contour3Display.PolarAxes.LastRadialAxisTextFontFile = ''
contour3Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''


# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------

anim = GetAnimationScene()
anim.UpdateAnimationUsingDataTimeSteps()
anim.GoToFirst()

import sys

nfr = anim.NumberOfFrames
nfr = len(fn)
for fr in range(nfr):
# save screenshot
  Log("Frame {:}/{:}".format(fr + 1, nfr))
  fn = 'aa.{:05d}.png'.format(fr + skipfirst)
  Log("Save to {:}".format(fn))
  SaveScreenshot(fn, renderView1, ImageResolution=[2000,2000])
  anim.GoToNext()

