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
    sys.stderr.write('''usage: {:} [traj_*.csv]
Plots isosurface.
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
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.OrientationAxesVisibility = 0
renderView1.CameraPosition =\
[1.673947199541654, 1.155534385955274, 3.017460244409912]
renderView1.CameraFocalPoint =\
[0.49676167169496643, -0.25002651930486053, -1.3449860440897603]
renderView1.CameraViewUp =\
[-0.0893981734829278, 0.9547941520217919, -0.2835067791833273]
renderView1.Background = [1.0]*3
renderView1.OSPRayMaterialLibrary = materialLibrary1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5

# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
traj = CSVReader(FileName=ff)

# list of all sources
vs = [traj]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
traj = ForceTime(traj)

# all ForceTime
vft = [traj]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# create a new 'Box'
box1 = Box()
box1.XLength = 2.0
box1.Center = [1.0, 0.5, 0.5]

table = TableToPoints(Input=traj)
table.XColumn = 'x'
table.YColumn = 'y'
table.ZColumn = 'z'

calc = Calculator(Input=table)
calc.AttributeType = 'Point Data'
calc.ResultArrayName = 'dz'
calc.Function = 'floor(coordsZ)'

warp = WarpByScalar(Input=calc)
warp.Scalars = ['POINTS', 'dz']
warp.Normal = [0.0, 0.0, -1.0]

glyph = Glyph(Input=warp, GlyphType='Sphere')
glyph.OrientationArray = ['POINTS', 'No orientation array']
glyph.ScaleArray = ['POINTS', 'r']
glyph.GlyphTransform = 'Transform2'
glyph.GlyphMode = 'All Points'
glyph.ScaleFactor = 1
glyph.GlyphType.Radius = 1.0
glyph.GlyphType.ThetaResolution = 9
glyph.GlyphType.PhiResolution = 9

thres = Threshold(Input=glyph)
thres.Scalars = ['POINTS', 'r']
thres.ThresholdRange = [0.0, 0.2]

box1Display = Show(box1, renderView1)
box1Display.Representation = 'Wireframe'
box1Display.AmbientColor = [0.0, 0.0, 0.0]
box1Display.ColorArrayName = [None, '']

thresDisplay = Show(thres, renderView1)
rLUT = GetColorTransferFunction('r')
rLUT.AutomaticRescaleRangeMode = 'Never'
rLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.030000000000000006, 0.865003, 0.865003, 0.865003, 0.06, 0.705882, 0.0156863, 0.14902]
rLUT.ScalarRangeInitialized = 1.0
rPWF = GetOpacityTransferFunction('r')
rPWF.Points = [0.0, 1.0, 0.5, 0.0, 0.06, 1.0, 0.5, 0.0]
rPWF.ScalarRangeInitialized = 1

thresDisplay.Representation = 'Surface'
thresDisplay.ColorArrayName = ['POINTS', 'r']
thresDisplay.LookupTable = rLUT
thresDisplay.ScalarOpacityFunction = rPWF


#####################################################
### END OF STATE FILE
#####################################################

for i in list(range(len(ss))):
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue

    with open(fn, "w") as f:
      f.write("")

    SetTime(i)

    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1)
