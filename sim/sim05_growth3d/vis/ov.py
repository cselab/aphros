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

# Create a new 'Light'
light1 = CreateLight()
light1.Intensity = 0.8
light1.Position = [4.0, 0.5, -8.0]
light1.FocalPoint = [4.0, 0.5, -6.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# create light
# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1200, 300]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [4.0, 0.498749986290932, 0.524999976158142]
renderView1.StereoType = 0
renderView1.CameraPosition = [7.05083065037415, 1.99328509007997, -15.0767923414364]
renderView1.CameraFocalPoint = [4.0, 0.498749986290932, 0.524999976158142]
renderView1.CameraViewUp = [-0.0010369827893544996, 0.9954619246947966, 0.09515503743694166]
renderView1.CameraParallelScale = 1.0564076736282666
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.196078431372549, 0.196078431372549, 0.196078431372549]
renderView1.EnableOSPRay = 1
renderView1.Shadows = 1
renderView1.AmbientSamples = 10
renderView1.SamplesPerPixel = 10
renderView1.AdditionalLights = light1
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

# create a new 'Box'
platey0 = Box()
platey0.XLength = 8.0
platey0.YLength = 0.05
platey0.Center = [4.0, -0.025, 0.5]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

def F(pre):
  l = GetFiles(pre)
  assert len(l) == nfr, "found %r files by '%r', expected %r" % (len(l), pre, nfr)
  return l[skipfirst:]


# create a new 'XDMF Reader'
fn = F("vf")
vf_ = XDMFReader(FileNames=fn)
vf_.CellArrayStatus = ['vf']
vf_.GridStatus = ['Grid_16802']

Log("Loaded {:} files from {:}".format(len(fn), base))

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=vf_)

# create a new 'Box'
electrode0 = Box()
electrode0.XLength = 0.05
electrode0.YLength = 0.4
electrode0.ZLength = 0.05
electrode0.Center = [0.525, -0.1975, 0.075]

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=vf_)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingDimensions = [512, 64, 67]
resampleToImage1.SamplingBounds = [0.0, 8.0, 0.0, 1.0, -0.05, 1.0]

# create a new 'Box'
electrode1 = Box()
electrode1.XLength = 0.05
electrode1.YLength = 0.4
electrode1.ZLength = 0.05
electrode1.Center = [0.525, 1.195, 0.075]

# create a new 'Append Geometry'
electrodes = AppendGeometry(Input=[electrode0, electrode1])

# create a new 'Contour'
contour1 = Contour(Input=resampleToImage1)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Clip'
oxygen = Clip(Input=contour1)
oxygen.ClipType = 'Plane'
oxygen.Scalars = ['POINTS', 'vtkValidPointMask']

# init the 'Plane' selected for 'ClipType'
oxygen.ClipType.Origin = [3.28865742683411, 0.5, 0.134546503424644]
oxygen.ClipType.Normal = [0.0, -1.0, 0.0]

# create a new 'Box'
platez1 = Box()
platez1.XLength = 8.0
platez1.YLength = 1.1
platez1.ZLength = 0.05
platez1.Center = [4.0, 0.5, 1.025]

# create a new 'Clip'
hydrogen = Clip(Input=contour1)
hydrogen.ClipType = 'Plane'
hydrogen.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
hydrogen.ClipType.Origin = [3.28865742683411, 0.5, 0.134546503424644]
hydrogen.ClipType.Normal = [0.0, 1.0, 0.0]

# create a new 'Box'
platez0 = Box()
platez0.XLength = 8.0
platez0.YLength = 1.1
platez0.ZLength = 0.05
platez0.Center = [4.0, 0.5, -0.025]

# create a new 'Box'
platey1 = Box()
platey1.XLength = 8.0
platey1.YLength = 0.05
platey1.Center = [4.0, 1.025, 0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from platey0
platey0Display = Show(platey0, renderView1)

# trace defaults for the display properties.
platey0Display.Representation = 'Surface'
platey0Display.ColorArrayName = [None, '']
platey0Display.DiffuseColor = [0.454901960784314, 0.454901960784314, 0.454901960784314]
platey0Display.OSPRayScaleArray = 'Normals'
platey0Display.OSPRayScaleFunction = 'PiecewiseFunction'
platey0Display.SelectOrientationVectors = 'None'
platey0Display.ScaleFactor = 0.8
platey0Display.SelectScaleArray = 'None'
platey0Display.GlyphType = 'Arrow'
platey0Display.GlyphTableIndexArray = 'None'
platey0Display.GaussianRadius = 0.00899999991059303
platey0Display.SetScaleArray = ['POINTS', 'Normals']
platey0Display.ScaleTransferFunction = 'PiecewiseFunction'
platey0Display.OpacityArray = ['POINTS', 'Normals']
platey0Display.OpacityTransferFunction = 'PiecewiseFunction'
platey0Display.DataAxesGrid = 'GridAxesRepresentation'
platey0Display.SelectionCellLabelFontFile = ''
platey0Display.SelectionPointLabelFontFile = ''
platey0Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
platey0Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
platey0Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
platey0Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.XTitleFontFile = ''
platey0Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.YTitleFontFile = ''
platey0Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.ZTitleFontFile = ''
platey0Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.XLabelFontFile = ''
platey0Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.YLabelFontFile = ''
platey0Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
platey0Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
platey0Display.PolarAxes.PolarAxisTitleFontFile = ''
platey0Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
platey0Display.PolarAxes.PolarAxisLabelFontFile = ''
platey0Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
platey0Display.PolarAxes.LastRadialAxisTextFontFile = ''
platey0Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
platey0Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from electrodes
electrodesDisplay = Show(electrodes, renderView1)

# trace defaults for the display properties.
electrodesDisplay.Representation = 'Surface'
electrodesDisplay.ColorArrayName = [None, '']
electrodesDisplay.DiffuseColor = [0.835294117647059, 0.556862745098039, 0.0]
electrodesDisplay.Diffuse = 0.63
electrodesDisplay.OSPRayScaleArray = 'Normals'
electrodesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
electrodesDisplay.SelectOrientationVectors = 'None'
electrodesDisplay.ScaleFactor = 0.8
electrodesDisplay.SelectScaleArray = 'None'
electrodesDisplay.GlyphType = 'Arrow'
electrodesDisplay.GlyphTableIndexArray = 'None'
electrodesDisplay.GaussianRadius = 0.00899999991059303
electrodesDisplay.SetScaleArray = ['POINTS', 'Normals']
electrodesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
electrodesDisplay.OpacityArray = ['POINTS', 'Normals']
electrodesDisplay.OpacityTransferFunction = 'PiecewiseFunction'
electrodesDisplay.DataAxesGrid = 'GridAxesRepresentation'
electrodesDisplay.SelectionCellLabelFontFile = ''
electrodesDisplay.SelectionPointLabelFontFile = ''
electrodesDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
electrodesDisplay.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
electrodesDisplay.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
electrodesDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.XTitleFontFile = ''
electrodesDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.YTitleFontFile = ''
electrodesDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.ZTitleFontFile = ''
electrodesDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.XLabelFontFile = ''
electrodesDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.YLabelFontFile = ''
electrodesDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
electrodesDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
electrodesDisplay.PolarAxes.PolarAxisTitleFontFile = ''
electrodesDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
electrodesDisplay.PolarAxes.PolarAxisLabelFontFile = ''
electrodesDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
electrodesDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
electrodesDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
electrodesDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from platey1
platey1Display = Show(platey1, renderView1)

# trace defaults for the display properties.
platey1Display.Representation = 'Surface'
platey1Display.AmbientColor = [0.0, 0.0, 0.0]
platey1Display.ColorArrayName = [None, '']
platey1Display.DiffuseColor = [0.454901960784314, 0.454901960784314, 0.454901960784314]
platey1Display.OSPRayScaleArray = 'Normals'
platey1Display.OSPRayScaleFunction = 'PiecewiseFunction'
platey1Display.SelectOrientationVectors = 'None'
platey1Display.ScaleFactor = 0.8
platey1Display.SelectScaleArray = 'None'
platey1Display.GlyphType = 'Arrow'
platey1Display.GlyphTableIndexArray = 'None'
platey1Display.GaussianRadius = 0.00899999991059303
platey1Display.SetScaleArray = ['POINTS', 'Normals']
platey1Display.ScaleTransferFunction = 'PiecewiseFunction'
platey1Display.OpacityArray = ['POINTS', 'Normals']
platey1Display.OpacityTransferFunction = 'PiecewiseFunction'
platey1Display.DataAxesGrid = 'GridAxesRepresentation'
platey1Display.SelectionCellLabelFontFile = ''
platey1Display.SelectionPointLabelFontFile = ''
platey1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
platey1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
platey1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
platey1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.XTitleFontFile = ''
platey1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.YTitleFontFile = ''
platey1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.ZTitleFontFile = ''
platey1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.XLabelFontFile = ''
platey1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.YLabelFontFile = ''
platey1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
platey1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
platey1Display.PolarAxes.PolarAxisTitleFontFile = ''
platey1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
platey1Display.PolarAxes.PolarAxisLabelFontFile = ''
platey1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
platey1Display.PolarAxes.LastRadialAxisTextFontFile = ''
platey1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
platey1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from platez1
platez1Display = Show(platez1, renderView1)

# trace defaults for the display properties.
platez1Display.Representation = 'Surface'
platez1Display.AmbientColor = [0.0, 0.0, 0.0]
platez1Display.ColorArrayName = [None, '']
platez1Display.DiffuseColor = [0.454901960784314, 0.454901960784314, 0.454901960784314]
platez1Display.OSPRayScaleArray = 'Normals'
platez1Display.OSPRayScaleFunction = 'PiecewiseFunction'
platez1Display.SelectOrientationVectors = 'None'
platez1Display.ScaleFactor = 0.8
platez1Display.SelectScaleArray = 'None'
platez1Display.GlyphType = 'Arrow'
platez1Display.GlyphTableIndexArray = 'None'
platez1Display.GaussianRadius = 0.00899999991059303
platez1Display.SetScaleArray = ['POINTS', 'Normals']
platez1Display.ScaleTransferFunction = 'PiecewiseFunction'
platez1Display.OpacityArray = ['POINTS', 'Normals']
platez1Display.OpacityTransferFunction = 'PiecewiseFunction'
platez1Display.DataAxesGrid = 'GridAxesRepresentation'
platez1Display.SelectionCellLabelFontFile = ''
platez1Display.SelectionPointLabelFontFile = ''
platez1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
platez1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
platez1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
platez1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.XTitleFontFile = ''
platez1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.YTitleFontFile = ''
platez1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.ZTitleFontFile = ''
platez1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.XLabelFontFile = ''
platez1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.YLabelFontFile = ''
platez1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
platez1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
platez1Display.PolarAxes.PolarAxisTitleFontFile = ''
platez1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
platez1Display.PolarAxes.PolarAxisLabelFontFile = ''
platez1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
platez1Display.PolarAxes.LastRadialAxisTextFontFile = ''
platez1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
platez1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from hydrogen
hydrogenDisplay = Show(hydrogen, renderView1)

# trace defaults for the display properties.
hydrogenDisplay.Representation = 'Surface'
hydrogenDisplay.AmbientColor = [0.0, 0.0, 0.0]
hydrogenDisplay.ColorArrayName = [None, '']
hydrogenDisplay.DiffuseColor = [0.541176470588235, 0.541176470588235, 0.541176470588235]
hydrogenDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
hydrogenDisplay.SelectOrientationVectors = 'None'
hydrogenDisplay.ScaleFactor = 0.549072283506393
hydrogenDisplay.SelectScaleArray = 'None'
hydrogenDisplay.GlyphType = 'Arrow'
hydrogenDisplay.GlyphTableIndexArray = 'None'
hydrogenDisplay.GaussianRadius = 0.0274536141753197
hydrogenDisplay.SetScaleArray = [None, '']
hydrogenDisplay.ScaleTransferFunction = 'PiecewiseFunction'
hydrogenDisplay.OpacityArray = [None, '']
hydrogenDisplay.OpacityTransferFunction = 'PiecewiseFunction'
hydrogenDisplay.DataAxesGrid = 'GridAxesRepresentation'
hydrogenDisplay.SelectionCellLabelFontFile = ''
hydrogenDisplay.SelectionPointLabelFontFile = ''
hydrogenDisplay.PolarAxes = 'PolarAxesRepresentation'
hydrogenDisplay.ScalarOpacityUnitDistance = 0.312750485185418

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
hydrogenDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.XTitleFontFile = ''
hydrogenDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.YTitleFontFile = ''
hydrogenDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.ZTitleFontFile = ''
hydrogenDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.XLabelFontFile = ''
hydrogenDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.YLabelFontFile = ''
hydrogenDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
hydrogenDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
hydrogenDisplay.PolarAxes.PolarAxisTitleFontFile = ''
hydrogenDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
hydrogenDisplay.PolarAxes.PolarAxisLabelFontFile = ''
hydrogenDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
hydrogenDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
hydrogenDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
hydrogenDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from oxygen
oxygenDisplay = Show(oxygen, renderView1)

# trace defaults for the display properties.
oxygenDisplay.Representation = 'Surface'
oxygenDisplay.AmbientColor = [0.0, 0.0, 0.0]
oxygenDisplay.ColorArrayName = [None, '']
oxygenDisplay.DiffuseColor = [0.541176470588235, 0.388235294117647, 0.392156862745098]
oxygenDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
oxygenDisplay.SelectOrientationVectors = 'None'
oxygenDisplay.ScaleFactor = 0.549072283506393
oxygenDisplay.SelectScaleArray = 'None'
oxygenDisplay.GlyphType = 'Arrow'
oxygenDisplay.GlyphTableIndexArray = 'None'
oxygenDisplay.GaussianRadius = 0.0274536141753197
oxygenDisplay.SetScaleArray = [None, '']
oxygenDisplay.ScaleTransferFunction = 'PiecewiseFunction'
oxygenDisplay.OpacityArray = [None, '']
oxygenDisplay.OpacityTransferFunction = 'PiecewiseFunction'
oxygenDisplay.DataAxesGrid = 'GridAxesRepresentation'
oxygenDisplay.SelectionCellLabelFontFile = ''
oxygenDisplay.SelectionPointLabelFontFile = ''
oxygenDisplay.PolarAxes = 'PolarAxesRepresentation'
oxygenDisplay.ScalarOpacityUnitDistance = 0.312750485185418

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
oxygenDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.XTitleFontFile = ''
oxygenDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.YTitleFontFile = ''
oxygenDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.ZTitleFontFile = ''
oxygenDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.XLabelFontFile = ''
oxygenDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.YLabelFontFile = ''
oxygenDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
oxygenDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
oxygenDisplay.PolarAxes.PolarAxisTitleFontFile = ''
oxygenDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
oxygenDisplay.PolarAxes.PolarAxisLabelFontFile = ''
oxygenDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
oxygenDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
oxygenDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
oxygenDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

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
  SaveScreenshot(fn, renderView1, ImageResolution=[4000,1000])
  anim.GoToNext()

