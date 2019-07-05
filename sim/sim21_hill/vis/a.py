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
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [1.0, 1.0, 1.0]
renderView1.CameraPosition = \
  [2.091994480544629, -5.193008444635638, 3.288843408789309]
renderView1.CameraFocalPoint = \
  [1.0, 1.0, 1.0]
renderView1.CameraViewUp = \
  [-0.059391174613884636, 0.33682408883346515, 0.9396926207859084]
renderView1.CameraParallelScale = 1.405915897498497
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.0, 1.0, 1.0]

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
clptom = CellDatatoPointData(Input=omm)

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

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
confvfDisplay.ScaleTransferFunction.Points = [-0.9999856352806091, 0.0, 0.5, 0.0, 0.9999890923500061, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
confvfDisplay.OpacityTransferFunction.Points = [-0.9999856352806091, 0.0, 0.5, 0.0, 0.9999890923500061, 1.0, 0.5, 0.0]

# create a new 'Resample To Image'
rsmp = ResampleToImage(Input=omm)
nx = 256
rsmp.SamplingDimensions = [nx + 1] * 3
dom = 2
rsmp.SamplingBounds = [0.0, dom, 0.0, dom, 0.0, dom]
calcomm = rsmp


# create a new 'Contour'
contour2 = Contour(Input=clptom)
contour2.ContourBy = ['POINTS', 'omm']
contour2.ComputeScalars = 1
contour2.Isosurfaces = np.linspace(0, 50, 10)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from appnd
appndDisplay = Show(omm, renderView1)
appndDisplay.Representation = 'Outline'
appndDisplay.ColorArrayName = ['CELLS', '']
appndDisplay.LineWidth = 2.0
appndDisplay.AmbientColor = [0.0, 0.0, 0.0]

# show data from contour2
contour2Display = Show(contour2, renderView1)
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.EnableOpacityMapping = 1
ommLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 15.000000000000004, 0.865003, 0.865003, 0.865003, 30.0, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0
contour2Display.Representation = 'Surface'
contour2Display.AmbientColor = [0.0, 0.0, 0.0]
contour2Display.ColorArrayName = ['POINTS', 'omm']
contour2Display.LookupTable = ommLUT
contour2Display.Opacity = 0.5
contour2Display.Ambient = 0.3

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [0.0, 0.17000000178813934, 0.5, 0.0, 30.0, 0.44999998807907104, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1

if not vort:
    Hide(contour2, renderView1)


#####################################################
### END OF STATE FILE
#####################################################

SetTime(1)
SaveScreenshot("tmp.png", renderView1)

for i in list(range(len(ss))):
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue

    SetTime(i)

    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1)

exit(0)
