# state file generated using paraview version 5.5.0-RC3

from paraview.simple import *

from glob import glob
import sys
import re
import math

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

def GetFiles(pre):
  global base
  p = base + "/{:}_*.csv".format(pre)
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

# create a new 'Box'
box1 = Box()
pi = math.pi
box1.XLength = 2. * pi
box1.YLength = 2. * pi
box1.ZLength = 2. * pi
box1.Center = [pi, pi, pi]

# create a new 'CSV Reader'
fn = F("traj")
csvreader = CSVReader(FileName=fn)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=csvreader)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Calculator'
calcr = Calculator(Input=tableToPoints1)
calcr.ResultArrayName = 'r'
calcr.Function = '(vf * 3 / 4 / 3.14)^(1/3)'

# create a new 'Calculator'
calcvel = Calculator(Input=tableToPoints1)
calcvel.ResultArrayName = 'vel'
calcvel.Function = 'vx*iHat+vy*jHat+vz*kHat'

# create a new 'Glyph'
glyph1 = Glyph(Input=calcvel,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'None']
glyph1.Vectors = ['POINTS', 'vel']
glyph1.ScaleMode = 'vector'
glyph1.ScaleFactor = 0.5
glyph1.GlyphMode = 'All Points'
glyph1.GlyphTransform = 'Transform2'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from box1
box1Display = Show(box1, renderView1)

# trace defaults for the display properties.
box1Display.Representation = 'Wireframe'
box1Display.ColorArrayName = [None, '']
box1Display.OSPRayScaleArray = 'Normals'
box1Display.OSPRayScaleFunction = 'PiecewiseFunction'
box1Display.SelectOrientationVectors = 'None'
box1Display.ScaleFactor = 0.628000020980835
box1Display.SelectScaleArray = 'None'
box1Display.GlyphType = 'Arrow'
box1Display.GlyphTableIndexArray = 'None'
box1Display.GaussianRadius = 0.03140000104904175
box1Display.SetScaleArray = ['POINTS', 'Normals']
box1Display.ScaleTransferFunction = 'PiecewiseFunction'
box1Display.OpacityArray = ['POINTS', 'Normals']
box1Display.OpacityTransferFunction = 'PiecewiseFunction'
box1Display.DataAxesGrid = 'GridAxesRepresentation'
box1Display.SelectionCellLabelFontFile = ''
box1Display.SelectionPointLabelFontFile = ''
box1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
box1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
box1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
box1Display.DataAxesGrid.XTitleFontFile = ''
box1Display.DataAxesGrid.YTitleFontFile = ''
box1Display.DataAxesGrid.ZTitleFontFile = ''
box1Display.DataAxesGrid.XLabelFontFile = ''
box1Display.DataAxesGrid.YLabelFontFile = ''
box1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
box1Display.PolarAxes.PolarAxisTitleFontFile = ''
box1Display.PolarAxes.PolarAxisLabelFontFile = ''
box1Display.PolarAxes.LastRadialAxisTextFontFile = ''
box1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from glyph1
glyph1Display = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.OSPRayScaleArray = 'GlyphVector'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 0.6983555465936662
glyph1Display.SelectScaleArray = 'GlyphVector'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'GlyphVector'
glyph1Display.GaussianRadius = 0.034917777329683305
glyph1Display.SetScaleArray = ['POINTS', 'GlyphVector']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'GlyphVector']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.SelectionCellLabelFontFile = ''
glyph1Display.SelectionPointLabelFontFile = ''
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-1.0501925945281982, 0.0, 0.5, 0.0, 1.1530451774597168, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-1.0501925945281982, 0.0, 0.5, 0.0, 1.1530451774597168, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleFontFile = ''
glyph1Display.DataAxesGrid.YTitleFontFile = ''
glyph1Display.DataAxesGrid.ZTitleFontFile = ''
glyph1Display.DataAxesGrid.XLabelFontFile = ''
glyph1Display.DataAxesGrid.YLabelFontFile = ''
glyph1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleFontFile = ''
glyph1Display.PolarAxes.PolarAxisLabelFontFile = ''
glyph1Display.PolarAxes.LastRadialAxisTextFontFile = ''
glyph1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from calcr
calcrDisplay = Show(calcr, renderView1)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.RGBPoints = [0.6525707835974837, 0.231373, 0.298039, 0.752941, 0.9967081099513457, 0.865003, 0.865003, 0.865003, 1.3408454363052076, 0.705882, 0.0156863, 0.14902]
pLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
calcrDisplay.Representation = 'Points'
calcrDisplay.ColorArrayName = ['POINTS', 'p']
calcrDisplay.LookupTable = pLUT
calcrDisplay.OSPRayUseScaleArray = 1
calcrDisplay.OSPRayScaleArray = 'r'
calcrDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
calcrDisplay.OSPRayMaterial = 'thin glass'
calcrDisplay.SelectOrientationVectors = 'None'
calcrDisplay.ScaleFactor = 0.6084381384713844
calcrDisplay.SelectScaleArray = 'r'
calcrDisplay.GlyphType = 'Arrow'
calcrDisplay.GlyphTableIndexArray = 'r'
calcrDisplay.GaussianRadius = 0.03042190692356922
calcrDisplay.SetScaleArray = ['POINTS', 'r']
calcrDisplay.ScaleTransferFunction = 'PiecewiseFunction'
calcrDisplay.OpacityArray = ['POINTS', 'r']
calcrDisplay.OpacityTransferFunction = 'PiecewiseFunction'
calcrDisplay.DataAxesGrid = 'GridAxesRepresentation'
calcrDisplay.SelectionCellLabelFontFile = ''
calcrDisplay.SelectionPointLabelFontFile = ''
calcrDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calcrDisplay.ScaleTransferFunction.Points = [0.09648003388216281, 0.0, 0.5, 0.0, 0.09715293552566363, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calcrDisplay.OpacityTransferFunction.Points = [0.09648003388216281, 0.0, 0.5, 0.0, 0.09715293552566363, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calcrDisplay.DataAxesGrid.XTitleFontFile = ''
calcrDisplay.DataAxesGrid.YTitleFontFile = ''
calcrDisplay.DataAxesGrid.ZTitleFontFile = ''
calcrDisplay.DataAxesGrid.XLabelFontFile = ''
calcrDisplay.DataAxesGrid.YLabelFontFile = ''
calcrDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calcrDisplay.PolarAxes.PolarAxisTitleFontFile = ''
calcrDisplay.PolarAxes.PolarAxisLabelFontFile = ''
calcrDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
calcrDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for pLUT in view renderView1
pLUTColorBar = GetScalarBar(pLUT, renderView1)
pLUTColorBar.Title = 'p'
pLUTColorBar.ComponentTitle = ''
pLUTColorBar.TitleFontFile = ''
pLUTColorBar.LabelFontFile = ''

# set color bar visibility
pLUTColorBar.Visibility = 1

# show color legend
calcrDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')
pPWF.Points = [0.6525707835974837, 0.0, 0.5, 0.0, 1.3408454363052076, 1.0, 0.5, 0.0]
pPWF.ScalarRangeInitialized = 1


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

