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

w = 1
wh = 2

def S(v):
  return [a * w for a in v]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [500, 1000 * wh]
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = S([2.785000000000005, 2.785000000000005, 11.140000000000063])
renderView1.CameraPosition = S([2.785000000000005, -74.27565179380603, 54.63098805708197])
renderView1.CameraFocalPoint = S([2.785000000000005, 2.785000000000005, 10.140000000000063])
renderView1.CameraViewUp = [0.0, 0.4999999999999998, 0.8660254037844388]
renderView1.CameraParallelScale = 5.134747553576296 * wh
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
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
confvfDisplay.DiffuseColor = [0.0, 1.0, 0.0]
confvfDisplay.AmbientColor = [0.0, 0.0, 0.0]

# create a new 'Append Attributes'
appendAttributes2 = AppendAttributes(Input=[omm])

# create a new 'Cell Data to Point Data'
#cellDatatoPointData2 = CellDatatoPointData(Input=appendAttributes2)

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=omm)
lim0, lim1 = GetBox(omm)
nz_nx = (lim1[2] - lim0[2]) / (lim1[0] - lim0[0])
nx = 128
resampleToImage1.SamplingDimensions = [nx + 1, nx + 1, int(nx * nz_nx + 0.5) + 1]
resampleToImage1.SamplingBounds = [l[i] for i in range(3) for l in [lim0,lim1]]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
resampleToImage1Display = Show(resampleToImage1, renderView1)

# get color transfer function/color map for 'omm'
om = 1.
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471,
                    om * 0.75, 0.865, 0.865, 0.865,
                    om * 1.5, 0.705882352941, 0.0156862745098, 0.149019607843]
ommLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [0.0, 0.0, 0.5, 0.0, om * 1.5, 0.28289473056793213, 0.5, 0.0]
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
resampleToImage1Display.SetScaleArray = ['POINTS', 'omm']
resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.OpacityArray = ['POINTS', 'omm']
resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleToImage1Display.SelectionCellLabelFontFile = ''
resampleToImage1Display.SelectionPointLabelFontFile = ''
resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
resampleToImage1Display.ScalarOpacityUnitDistance = 0.11630404359472776
resampleToImage1Display.ScalarOpacityFunction = ommPWF
resampleToImage1Display.IsosurfaceValues = [0.5000000000000001]
resampleToImage1Display.Slice = 31
resampleToImage1Display.Shade = 1

resampleToImage1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
resampleToImage1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
resampleToImage1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

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
