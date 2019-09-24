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
renderView1.ViewSize = [1024, 1024]
renderView1.OrientationAxesVisibility = 0

renderView1.CameraPosition = \
  [1.8523389696485202, 1.3297883193504907, 2.801482605769909]
renderView1.CameraFocalPoint = \
  [0.20376001199494642, 0.22701851625700403, 0.10662850870392106]
renderView1.CameraViewUp = \
  [-0.1701681903046653, 0.9441282598262686, -0.28224920904381173]
renderView1.CameraParallelScale = \
  0.866025403784

if 1: # extent 2
  renderView1.CameraPosition = \
    [0.5, 3.3660254037844375, 4.598076211353318]
  renderView1.CameraFocalPoint = \
    [0.5, 1.0, 0.5]
  renderView1.CameraViewUp = \
    [0.0, 0.8660254037844387, -0.4999999999999997]
  renderView1.CameraParallelScale = \
    1.15633239565
  renderView1.CameraParallelProjection = 1

renderView1.Shadows = 0
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.EnableOSPRay = 0
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 5
renderView1.OSPRayMaterialLibrary = materialLibrary1

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

if 0:
  clip1 = Clip(Input=confvf)
  clip1.ClipType = 'Plane'
  clip1.ClipType.Origin = [0.5, 0.5, 0.5]
  clip1.ClipType.Normal = [1.0, 0.0, 1.0]
  confvf = clip1

confvfDisplay = Show(confvf, renderView1)
confvfDisplay.Representation = 'Surface'
confvfDisplay.ColorArrayName = [None, '']
confvfDisplay.Opacity = 0.5


calcommDisplay = Show(omm, renderView1)
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [5.0, 0.231373, 0.298039, 0.752941, 52.50000000000001, 0.865003, 0.865003, 0.865003, 100.0, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [5.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1


calcommDisplay.Representation = 'Volume'
calcommDisplay.AmbientColor = [0.0, 0.0, 0.0]
calcommDisplay.ColorArrayName = ['CELLS', 'omm']
calcommDisplay.LookupTable = ommLUT
calcommDisplay.ScalarOpacityUnitDistance = 0.04251092259923938
calcommDisplay.ScalarOpacityFunction = ommPWF


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
