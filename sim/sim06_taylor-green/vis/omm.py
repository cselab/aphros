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


ospray = 1
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
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [3.14159274101257, 3.14159274101257, 3.14159274101257]
renderView1.StereoType = 0
renderView1.CameraPosition = [8.535405181306203, -11.84729921683567, 8.958582336091343]
renderView1.CameraFocalPoint = [1.7784384184717228, 6.717314388057494, 1.7679684977243724]
renderView1.CameraViewUp = [-0.11697777844051055, 0.3213938048432693, 0.9396926207859086]
renderView1.CameraParallelScale = 4.596149164737374
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 5
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

clpt = CellDatatoPointData(Input=vf)

# create a new 'Contour'
confvf = Contour(Input=clpt)
confvf.ContourBy = ['POINTS', 'vf']
confvf.Isosurfaces = [0.5]
confvf.PointMergeMethod = 'Uniform Binning'

# show data from confvf
confvfDisplay = Show(confvf, renderView1)

# trace defaults for the display properties.
confvfDisplay.Representation = 'Surface'
confvfDisplay.ColorArrayName = [None, '']
confvfDisplay.LineWidth = 3.0
confvfDisplay.OSPRayScaleArray = 'Normals'
confvfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
confvfDisplay.SelectOrientationVectors = 'None'
confvfDisplay.ScaleFactor = 0.1
confvfDisplay.SelectScaleArray = 'None'
confvfDisplay.GlyphType = 'Arrow'
confvfDisplay.GlyphTableIndexArray = 'None'
confvfDisplay.GaussianRadius = 0.005
confvfDisplay.SetScaleArray = ['POINTS', 'Normals']
confvfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
confvfDisplay.OpacityArray = ['POINTS', 'Normals']
confvfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
confvfDisplay.DataAxesGrid = 'GridAxesRepresentation'
confvfDisplay.SelectionCellLabelFontFile = ''
confvfDisplay.SelectionPointLabelFontFile = ''
confvfDisplay.PolarAxes = 'PolarAxesRepresentation'
confvfDisplay.DiffuseColor = [1.0, 0.0, 0.0]
confvfDisplay.AmbientColor = [0.0, 0.0, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
confvfDisplay.ScaleTransferFunction.Points = [-0.9999856352806091, 0.0, 0.5, 0.0, 0.9999890923500061, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
confvfDisplay.OpacityTransferFunction.Points = [-0.9999856352806091, 0.0, 0.5, 0.0, 0.9999890923500061, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
confvfDisplay.DataAxesGrid.XTitleFontFile = ''
confvfDisplay.DataAxesGrid.YTitleFontFile = ''
confvfDisplay.DataAxesGrid.ZTitleFontFile = ''
confvfDisplay.DataAxesGrid.XLabelFontFile = ''
confvfDisplay.DataAxesGrid.YLabelFontFile = ''
confvfDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
confvfDisplay.PolarAxes.PolarAxisTitleFontFile = ''
confvfDisplay.PolarAxes.PolarAxisLabelFontFile = ''
confvfDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
confvfDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# create a new 'Resample To Image'
rsmp = ResampleToImage(Input=omm)
nx = 256
rsmp.SamplingDimensions = [nx + 1] * 3
pi = 3.141592
dom = 2 * pi
rsmp.SamplingBounds = [0.0, dom, 0.0, dom, 0.0, dom]

calcomm = rsmp

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from appnd
appndDisplay = Show(omm, renderView1)

# trace defaults for the display properties.
appndDisplay.Representation = 'Outline'
appndDisplay.ColorArrayName = ['CELLS', '']
appndDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
appndDisplay.SelectOrientationVectors = 'None'
appndDisplay.ScaleFactor = 0.30000000000000004
appndDisplay.SelectScaleArray = 'vf'
appndDisplay.GlyphType = 'Arrow'
appndDisplay.GlyphTableIndexArray = 'vf'
appndDisplay.GaussianRadius = 0.015
appndDisplay.SetScaleArray = [None, '']
appndDisplay.ScaleTransferFunction = 'PiecewiseFunction'
appndDisplay.OpacityArray = [None, '']
appndDisplay.OpacityTransferFunction = 'PiecewiseFunction'
appndDisplay.DataAxesGrid = 'GridAxesRepresentation'
appndDisplay.SelectionCellLabelFontFile = ''
appndDisplay.SelectionPointLabelFontFile = ''
appndDisplay.PolarAxes = 'PolarAxesRepresentation'
appndDisplay.AmbientColor = [0.0, 0.0, 0.0]


# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
appndDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
appndDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
appndDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
appndDisplay.DataAxesGrid.XTitleFontFile = ''
appndDisplay.DataAxesGrid.YTitleFontFile = ''
appndDisplay.DataAxesGrid.ZTitleFontFile = ''
appndDisplay.DataAxesGrid.XLabelFontFile = ''
appndDisplay.DataAxesGrid.YLabelFontFile = ''
appndDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
appndDisplay.PolarAxes.PolarAxisTitleFontFile = ''
appndDisplay.PolarAxes.PolarAxisLabelFontFile = ''
appndDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
appndDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from calcomm
calcommDisplay = Show(calcomm, renderView1)

if not vort:
    Hide(calcoom, renderView1)

k = 0.5
# get color transfer function/color map for 'omm'
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [
    1.0 * k, 0.231373, 0.298039, 0.752941,
    5.5 * k, 0.865003, 0.865003, 0.865003,
    10.0 * k, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [
    1.0 * k, 0.0, 0.5, 0.0,
    10.0 * k, 0.5073529481887817, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1


# trace defaults for the display properties.
calcommDisplay.Representation = 'Volume'
calcommDisplay.AmbientColor = [0.0, 0.0, 0.0]
calcommDisplay.ColorArrayName = ['POINTS', 'omm']
calcommDisplay.LookupTable = ommLUT
calcommDisplay.OSPRayScaleArray = 'omm'
calcommDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
calcommDisplay.SelectOrientationVectors = 'vel'
calcommDisplay.ScaleFactor = 0.6283185307179574
calcommDisplay.SelectScaleArray = 'omm'
calcommDisplay.GlyphType = 'Arrow'
calcommDisplay.GlyphTableIndexArray = 'omm'
calcommDisplay.GaussianRadius = 0.03141592653589787
calcommDisplay.SetScaleArray = ['POINTS', 'omm']
calcommDisplay.ScaleTransferFunction = 'PiecewiseFunction'
calcommDisplay.OpacityArray = ['POINTS', 'omm']
calcommDisplay.OpacityTransferFunction = 'PiecewiseFunction'
calcommDisplay.DataAxesGrid = 'GridAxesRepresentation'
calcommDisplay.SelectionCellLabelFontFile = ''
calcommDisplay.SelectionPointLabelFontFile = ''
calcommDisplay.PolarAxes = 'PolarAxesRepresentation'
calcommDisplay.ScalarOpacityUnitDistance = 0.1700436903969576
calcommDisplay.ScalarOpacityFunction = ommPWF
calcommDisplay.IsosurfaceValues = [7.225581137787295]
calcommDisplay.Slice = 32

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calcommDisplay.ScaleTransferFunction.Points = [0.009214281748931607, 0.0, 0.5, 0.0, 13.063280207727544, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calcommDisplay.OpacityTransferFunction.Points = [0.009214281748931607, 0.0, 0.5, 0.0, 13.063280207727544, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

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
