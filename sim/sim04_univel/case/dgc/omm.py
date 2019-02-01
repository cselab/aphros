#!/usr/bin/env pvbatch

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

# Sets time of datasets to step i
def SetTime(i):
    global vft, vt
    for j in range(len(vft)):
        s = vft[j]
        s.ForcedTime = vt[j][i]
        s.UpdatePipeline()

# Returns bounding box of object o
def GetBox(o):
    o.UpdatePipeline()
    di = o.GetDataInformation()
    lim = di.DataInformation.GetBounds()
    lim0 = np.array(lim[::2])
    lim1 = np.array(lim[1::2])
    return lim0, lim1


ospray = 0
vort = 1

av = sys.argv
if len(av) < 2 or av[1] == '-h':
    sys.stderr.write('''usage: {:} [vf_*.xmf]
Plots bubbles.
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

# vf input
ff = natsorted(av[1:])
# vf basename
ffb = list(map(os.path.basename, ff))
# vf dirname
ffd = list(map(os.path.dirname, ff))
# steps
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]
# omm input
ffomm = ["omm_{:04d}.xmf".format(s) for s in ss]
# append dirname
for i in range(len(ss)):
  ffomm[i] = os.path.join(ffd[i], ffomm[i])


# output pattern (:0 substituted by frame number)
bo = "a_{:}.png"

#####################################################
### BEGIN OF STATE FILE
#####################################################

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# RENDER
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000, 1000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 1.0, 0.5]
renderView1.StereoType = 0

def T(a):
  return [a[0], a[2], a[1]]

renderView1.CameraPosition = T([3.5378747228904905, -16.79505969829756, 6.888169545748997])
renderView1.CameraFocalPoint = T([0.10727278958569013, 2.660850676399359, -0.3024442926179669])
renderView1.CameraViewUp = T([-0.05939117461388435, 0.3368240888334645, 0.9396926207859086])
renderView1.CameraParallelScale = 0.6823124917434552
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.EnableOSPRay = ospray
renderView1.AmbientSamples = 10
renderView1.SamplesPerPixel = 10
renderView1.Shadows = 1
renderView1.OSPRayMaterialLibrary = materialLibrary1

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
vf = XDMFReader(FileNames=ff)
vf.CellArrayStatus = ['vf']
vf.GridStatus = ['Grid_0']


omm = XDMFReader(FileNames=ffomm)
omm.CellArrayStatus = ['omm']
omm.GridStatus = ['Grid_1']

# list of all sources
vs = [vf, omm]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
vf = ForceTime(vf)
omm = ForceTime(omm)

# all ForceTime
vft = [vf, omm]

SetTime(0)

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# create a new 'Append Attributes'
appendAttributes2 = AppendAttributes(Input=[omm, vf])

# create a new 'Cell Data to Point Data'
cellDatatoPointData2 = CellDatatoPointData(Input=appendAttributes2)

# create a new 'Contour'
contour2 = Contour(Input=cellDatatoPointData2)
contour2.ContourBy = ['POINTS', 'vf']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=cellDatatoPointData2)
s = omm
nx0,nx1,ny0,ny1,nz0,nz1 = omm.GetDataInformation().GetExtent()
bb = omm.GetDataInformation().GetBounds()
nx = nx1 - nx0 + 1
ny = ny1 - ny0 + 1
nz = nz1 - nz0 + 1
print("Size: {:} {:} {:}".format(nx, ny, nz))
print("Bounds: {:}".format(bb))
resampleToImage1.SamplingDimensions = [nx, ny, nz]
resampleToImage1.SamplingBounds = bb

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from contour2
contour2Display = Show(contour2, renderView1)

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = [None, '']
contour2Display.Opacity = 0.75
contour2Display.OSPRayScaleArray = 'Normals'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.1
contour2Display.SelectScaleArray = 'None'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'None'
contour2Display.GaussianRadius = 0.005
contour2Display.SetScaleArray = ['POINTS', 'Normals']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'Normals']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.SelectionCellLabelFontFile = ''
contour2Display.SelectionPointLabelFontFile = ''
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour2Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [-0.9978506565093994, 1.0, 0.5, 0.0, 0.9989058971405029, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [-0.9978506565093994, 1.0, 0.5, 0.0, 0.9989058971405029, 1.0, 0.5, 0.0]

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

# show data from resampleToImage1
resampleToImage1Display = Show(resampleToImage1, renderView1)

# get color transfer function/color map for 'omm'
omk = 3.
ommLUT = GetColorTransferFunction('omm')
ommLUT.RGBPoints = [0, 0.231373, 0.298039, 0.752941, omk * 2., 0.865003, 0.865003, 0.865003, omk * 4, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [0, 0.0, 0.5, 0.0, omk * 0.2, 0.0, 0.5, 0.0, omk * 4, 0.5, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
resampleToImage1Display.Representation = 'Volume'
resampleToImage1Display.ColorArrayName = ['POINTS', 'omm']
resampleToImage1Display.LookupTable = ommLUT
resampleToImage1Display.OSPRayScaleArray = 'vf'
resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
resampleToImage1Display.SelectOrientationVectors = 'None'
resampleToImage1Display.ScaleFactor = 0.2
resampleToImage1Display.SelectScaleArray = 'vf'
resampleToImage1Display.GlyphType = 'Arrow'
resampleToImage1Display.GlyphTableIndexArray = 'vf'
resampleToImage1Display.GaussianRadius = 0.01
resampleToImage1Display.SetScaleArray = ['POINTS', 'vf']
resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.OpacityArray = ['POINTS', 'vf']
resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleToImage1Display.SelectionCellLabelFontFile = ''
resampleToImage1Display.SelectionPointLabelFontFile = ''
resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
resampleToImage1Display.ScalarOpacityUnitDistance = 0.030377520269369778
resampleToImage1Display.ScalarOpacityFunction = ommPWF
resampleToImage1Display.IsosurfaceValues = [0.5000000000000001]
resampleToImage1Display.Slice = 31

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
resampleToImage1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
resampleToImage1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
resampleToImage1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

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

# show data from appendAttributes2
appendAttributes2Display = Show(appendAttributes2, renderView1)

# trace defaults for the display properties.
appendAttributes2Display.Representation = 'Outline'
appendAttributes2Display.ColorArrayName = ['CELLS', '']
appendAttributes2Display.LineWidth = 3.0
appendAttributes2Display.RenderLinesAsTubes = 1
appendAttributes2Display.OSPRayScaleFunction = 'PiecewiseFunction'
appendAttributes2Display.SelectOrientationVectors = 'None'
appendAttributes2Display.ScaleFactor = 0.2
appendAttributes2Display.SelectScaleArray = 'omm'
appendAttributes2Display.GlyphType = 'Arrow'
appendAttributes2Display.GlyphTableIndexArray = 'omm'
appendAttributes2Display.GaussianRadius = 0.01
appendAttributes2Display.SetScaleArray = [None, '']
appendAttributes2Display.ScaleTransferFunction = 'PiecewiseFunction'
appendAttributes2Display.OpacityArray = [None, '']
appendAttributes2Display.OpacityTransferFunction = 'PiecewiseFunction'
appendAttributes2Display.DataAxesGrid = 'GridAxesRepresentation'
appendAttributes2Display.SelectionCellLabelFontFile = ''
appendAttributes2Display.SelectionPointLabelFontFile = ''
appendAttributes2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
appendAttributes2Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
appendAttributes2Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
appendAttributes2Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
appendAttributes2Display.DataAxesGrid.XTitleFontFile = ''
appendAttributes2Display.DataAxesGrid.YTitleFontFile = ''
appendAttributes2Display.DataAxesGrid.ZTitleFontFile = ''
appendAttributes2Display.DataAxesGrid.XLabelFontFile = ''
appendAttributes2Display.DataAxesGrid.YLabelFontFile = ''
appendAttributes2Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
appendAttributes2Display.PolarAxes.PolarAxisTitleFontFile = ''
appendAttributes2Display.PolarAxes.PolarAxisLabelFontFile = ''
appendAttributes2Display.PolarAxes.LastRadialAxisTextFontFile = ''
appendAttributes2Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''


#####################################################
### END OF STATE FILE
#####################################################

for i in list(range(len(ss))):
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue

    SetTime(i)

    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1)

exit(0)
