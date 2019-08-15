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
Plots bubbles with vorticity.
Current folder:
omm_*.xmf
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

# input
ff = natsorted(av[1:])
# basename
ffb = list(map(os.path.basename, ff))
# dirname
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

# width
w=5.57

def S(v):
  return [a * w for a in v]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [700, 1400]
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = S([0.5, 0.5, 1.0])
renderView1.CameraPosition = S(
    [8.023735598478023, 3.3063072646186735, 3.641632333070597])
renderView1.CameraFocalPoint = S(
    [0.5, 0.5, 2.0])
renderView1.CameraViewUp = S(
    [-0.18725081813571506, -0.07110396830462755, 0.9797353503874605])
renderView1.CameraParallelScale = 1.9159375215828818
#renderView1.CameraParallelProjection = 1
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

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# create a new 'Contour'
clpt = CellDatatoPointData(Input=vf)
confvf = Contour(Input=clpt)
confvf.ContourBy = ['POINTS', 'vf']
confvf.Isosurfaces = [0.5]
confvf.PointMergeMethod = 'Uniform Binning'
confvfDisplay = Show(confvf, renderView1)
confvfDisplay.Representation = 'Surface'
confvfDisplay.ColorArrayName = [None, '']

# create a new 'Append Attributes'
appendAttributes2 = AppendAttributes(Input=[omm])

# create a new 'Cell Data to Point Data'
cellDatatoPointData2 = CellDatatoPointData(Input=appendAttributes2)

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=cellDatatoPointData2)
resampleToImage1.SamplingDimensions = [129, 129, 257]
resampleToImage1.SamplingBounds = [0.0, w, 0.0, w, 0.0, w * 2]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
resampleToImage1Display = Show(resampleToImage1, renderView1)

# get color transfer function/color map for 'omm'
ommLUT = GetColorTransferFunction('omm')
om = 1.
ommLUT.RGBPoints = [0., 0.231373, 0.298039, 0.752941,
                    om, 0.865003, 0.865003, 0.865003,
                    om * 2, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [0., 0, 0.5, 0.,
                 om, 0.25, 0.5, 0.,
                 om * 2, 0.5, 0.5, 0.]
ommPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
resampleToImage1Display.Representation = 'Volume'
resampleToImage1Display.ColorArrayName = ['POINTS', 'omm']
resampleToImage1Display.LookupTable = ommLUT
resampleToImage1Display.SelectOrientationVectors = 'None'
resampleToImage1Display.ScaleFactor = 0.2
resampleToImage1Display.GlyphType = 'Arrow'
resampleToImage1Display.GaussianRadius = 0.01
resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
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

'''
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
'''


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
