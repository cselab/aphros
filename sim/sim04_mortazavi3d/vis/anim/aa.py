# state file generated using paraview version 5.5.0-RC3

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.0-RC3

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

log = open("log", 'w')

from glob import glob
import sys
import re

# natural sort
def natkey(s, _nsre=re.compile('([0-9]+)')):
      return [int(text) if text.isdigit() else text.lower()
                      for text in re.split(_nsre, s)]

def natsorted(v):
  return sorted(v, key=natkey)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [8.0, 0.5, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [32.2629821062112, 5.14754209840159, -18.2770381080199]
renderView1.CameraFocalPoint = [8.0, 0.5, 0.499999999999999]
renderView1.CameraViewUp = [-0.11853115278243698, 0.9887200778001208, 0.09155858001842139]
renderView1.CameraParallelScale = 3.25
renderView1.CameraParallelProjection = 1
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
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
base = "/users/karnakov/s/scratch/ch/bub/sim04_p2048_lx16_t20/"
#base = "/users/karnakov/s/scratch/ch/bub/sim04_p32/"

# create a new 'XDMF Reader'
fn = natsorted(glob(base + "vx_*.xmf"))
vx_xmf = XDMFReader(FileNames=fn)
vx_xmf.CellArrayStatus = ['vx']
vx_xmf.GridStatus = ['Grid_16805']

# create a new 'XDMF Reader'
fn = natsorted(glob(base + "vy_*.xmf"))
vy_xmf = XDMFReader(FileNames=fn)
vy_xmf.CellArrayStatus = ['vy']
vy_xmf.GridStatus = ['Grid_16808']

fn = natsorted(glob(base + "vz_*.xmf"))
vz_xmf = XDMFReader(FileNames=fn)
vz_xmf.CellArrayStatus = ['vz']
vz_xmf.GridStatus = ['Grid_16811']

# create a new 'XDMF Reader'
fn = natsorted(glob(base + "vf_*.xmf"))
vf_xmf = XDMFReader(FileNames=fn)
vf_xmf.CellArrayStatus = ['vf']
vf_xmf.GridStatus = ['Grid_16802']

log.write("Loaded {:} files from {:}\n".format(len(fn), base))



# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(Input=[vf_xmf, vx_xmf, vy_xmf, vz_xmf])

# create a new 'Calculator'
calculator1 = Calculator(Input=appendAttributes1)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'vel'
calculator1.Function = 'vx*iHat+vy*jHat+vz*kHat'

# create a new 'Gradient Of Unstructured DataSet'
gradientOfUnstructuredDataSet1 = GradientOfUnstructuredDataSet(Input=calculator1)
gradientOfUnstructuredDataSet1.ScalarArray = ['CELLS', 'vel']
gradientOfUnstructuredDataSet1.ComputeGradient = 0
gradientOfUnstructuredDataSet1.ComputeVorticity = 1

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=gradientOfUnstructuredDataSet1)
resampleToImage1.SamplingDimensions = [2048, 128, 128]
resampleToImage1.SamplingBounds = [0.0, 4.0, 0.0, 1.0, 0.0, 1.0]

# create a new 'Contour'
contour1 = Contour(Input=resampleToImage1)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Slice'
slice1 = Slice(Input=resampleToImage1)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [-4.0, 4.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [8.0, 0.5, 0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from gradientOfUnstructuredDataSet1
gradientOfUnstructuredDataSet1Display = Show(gradientOfUnstructuredDataSet1, renderView1)

# trace defaults for the display properties.
gradientOfUnstructuredDataSet1Display.Representation = 'Outline'
gradientOfUnstructuredDataSet1Display.ColorArrayName = ['CELLS', '']
gradientOfUnstructuredDataSet1Display.PointSize = 1.0
gradientOfUnstructuredDataSet1Display.LineWidth = 3.0
gradientOfUnstructuredDataSet1Display.RenderLinesAsTubes = 1
gradientOfUnstructuredDataSet1Display.OSPRayScaleFunction = 'PiecewiseFunction'
gradientOfUnstructuredDataSet1Display.SelectOrientationVectors = 'vel'
gradientOfUnstructuredDataSet1Display.ScaleFactor = 0.4
gradientOfUnstructuredDataSet1Display.SelectScaleArray = 'vf'
gradientOfUnstructuredDataSet1Display.GlyphType = 'Arrow'
gradientOfUnstructuredDataSet1Display.GlyphTableIndexArray = 'vf'
gradientOfUnstructuredDataSet1Display.GaussianRadius = 0.02
gradientOfUnstructuredDataSet1Display.SetScaleArray = [None, '']
gradientOfUnstructuredDataSet1Display.ScaleTransferFunction = 'PiecewiseFunction'
gradientOfUnstructuredDataSet1Display.OpacityArray = [None, '']
gradientOfUnstructuredDataSet1Display.OpacityTransferFunction = 'PiecewiseFunction'
gradientOfUnstructuredDataSet1Display.DataAxesGrid = 'GridAxesRepresentation'
gradientOfUnstructuredDataSet1Display.SelectionCellLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.SelectionPointLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
gradientOfUnstructuredDataSet1Display.DataAxesGrid.XTitleFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.YTitleFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.ZTitleFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.XLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.YLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
gradientOfUnstructuredDataSet1Display.PolarAxes.PolarAxisTitleFontFile = ''
gradientOfUnstructuredDataSet1Display.PolarAxes.PolarAxisLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.PolarAxes.LastRadialAxisTextFontFile = ''
gradientOfUnstructuredDataSet1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from resampleToImage1
resampleToImage1Display = Show(resampleToImage1, renderView1)

# get color transfer function/color map for 'Vorticity'
vorticityLUT = GetColorTransferFunction('Vorticity')
vorticityLUT.AutomaticRescaleRangeMode = 'Never'
vorticityLUT.RGBPoints = [0.0, 0.0416667, 0.0, 0.0, 1.90476, 0.208333, 0.0, 0.0, 3.80952, 0.375, 0.0, 0.0, 5.71428, 0.541667, 0.0, 0.0, 7.619055, 0.708333, 0.0, 0.0, 9.52381500000003, 0.854137, 0.0, 0.0, 11.428575, 0.937488, 0.039062, 0.0, 13.333335, 1.0, 0.208333, 0.0, 15.238095, 1.0, 0.375, 0.0, 17.142855, 1.0, 0.541667, 0.0, 19.047615, 1.0, 0.708333, 0.0, 20.952375, 1.0, 0.858805, 0.03125, 22.85715, 1.0, 0.947392, 0.15625, 24.76191, 1.0, 1.0, 0.3125, 26.66667, 1.0, 1.0, 0.5625, 28.57143, 1.0, 1.0, 0.8125, 30.0, 1.0, 1.0, 1.0]
vorticityLUT.ColorSpace = 'Lab'
vorticityLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Vorticity'
vorticityPWF = GetOpacityTransferFunction('Vorticity')
vorticityPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.74390250444412, 0.0, 0.5, 0.0, 8.78048801422118, 0.394736856222153, 0.5, 0.0, 30.0, 0.717105269432068, 0.5, 0.0]
vorticityPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
resampleToImage1Display.Representation = 'Volume'
resampleToImage1Display.ColorArrayName = ['POINTS', 'Vorticity']
resampleToImage1Display.LookupTable = vorticityLUT
resampleToImage1Display.OSPRayScaleArray = 'Vorticity'
resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
resampleToImage1Display.SelectOrientationVectors = 'None'
resampleToImage1Display.ScaleFactor = 0.4
resampleToImage1Display.SelectScaleArray = 'None'
resampleToImage1Display.GlyphType = 'Arrow'
resampleToImage1Display.GlyphTableIndexArray = 'None'
resampleToImage1Display.GaussianRadius = 0.02
resampleToImage1Display.SetScaleArray = ['POINTS', 'Vorticity']
resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.OpacityArray = ['POINTS', 'Vorticity']
resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleToImage1Display.SelectionCellLabelFontFile = ''
resampleToImage1Display.SelectionPointLabelFontFile = ''
resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
resampleToImage1Display.ScalarOpacityUnitDistance = 0.0210035870888237
resampleToImage1Display.ScalarOpacityFunction = vorticityPWF
resampleToImage1Display.Slice = 63

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
resampleToImage1Display.ScaleTransferFunction.Points = [-222.232727663729, 0.0, 0.5, 0.0, 223.847467804814, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
resampleToImage1Display.OpacityTransferFunction.Points = [-222.232727663729, 0.0, 0.5, 0.0, 223.847467804814, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
resampleToImage1Display.DataAxesGrid.XTitleFontFile = ''
resampleToImage1Display.DataAxesGrid.YTitleFontFile = ''
resampleToImage1Display.DataAxesGrid.ZTitleFontFile = ''
resampleToImage1Display.DataAxesGrid.XLabelFontFile = ''
resampleToImage1Display.DataAxesGrid.YLabelFontFile = ''
resampleToImage1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
resampleToImage1Display.PolarAxes.PolarAxisTitleFontFile = ''
resampleToImage1Display.PolarAxes.PolarAxisLabelFontFile = ''
resampleToImage1Display.PolarAxes.LastRadialAxisTextFontFile = ''
resampleToImage1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']
contour1Display.OSPRayScaleArray = 'Normals'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.307632088661194
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.GaussianRadius = 0.0153816044330597
contour1Display.SetScaleArray = ['POINTS', 'Normals']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'Normals']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.SelectionCellLabelFontFile = ''
contour1Display.SelectionPointLabelFontFile = ''
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-0.999942541122437, 0.0, 0.5, 0.0, 0.999990284442902, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-0.999942541122437, 0.0, 0.5, 0.0, 0.999990284442902, 1.0, 0.5, 0.0]

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

# show data from slice1
slice1Display = Show(slice1, renderView1)

# get color transfer function/color map for 'vel'
velLUT = GetColorTransferFunction('vel')
velLUT.RGBPoints = [0.000424215082919017, 0.231373, 0.298039, 0.752941, 0.277055692278943, 0.865003, 0.865003, 0.865003, 0.553687169474967, 0.705882, 0.0156863, 0.14902]
velLUT.ScalarRangeInitialized = 1.0
velLUT.VectorMode = 'Component'

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'vel']
slice1Display.LookupTable = velLUT
slice1Display.Ambient = 0.36
slice1Display.OSPRayScaleArray = 'Vorticity'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.8
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.04
slice1Display.SetScaleArray = ['POINTS', 'Vorticity']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'Vorticity']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.SelectionCellLabelFontFile = ''
slice1Display.SelectionPointLabelFontFile = ''
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [-0.762280949621196, 0.0, 0.5, 0.0, 0.722671741484164, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [-0.762280949621196, 0.0, 0.5, 0.0, 0.722671741484164, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
slice1Display.DataAxesGrid.XTitleFontFile = ''
slice1Display.DataAxesGrid.YTitleFontFile = ''
slice1Display.DataAxesGrid.ZTitleFontFile = ''
slice1Display.DataAxesGrid.XLabelFontFile = ''
slice1Display.DataAxesGrid.YLabelFontFile = ''
slice1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
slice1Display.PolarAxes.PolarAxisTitleFontFile = ''
slice1Display.PolarAxes.PolarAxisLabelFontFile = ''
slice1Display.PolarAxes.LastRadialAxisTextFontFile = ''
slice1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'vel'
velPWF = GetOpacityTransferFunction('vel')
velPWF.Points = [0.000424215082919017, 0.0, 0.5, 0.0, 0.553687169474967, 1.0, 0.5, 0.0]
velPWF.ScalarRangeInitialized = 1

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
  log.write("Frame {:}/{:}\n".format(fr + 1, nfr))
  log.flush()
  SaveScreenshot('aa.{:05d}.png'.format(fr), renderView1, ImageResolution=[1920, 1080])
  anim.GoToNext()

