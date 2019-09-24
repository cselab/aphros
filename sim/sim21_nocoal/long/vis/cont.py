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
renderView1.ViewSize = [1400, 700]
renderView1.OrientationAxesVisibility = 0

renderView1.CameraPosition = \
  [0.8828178821634488, -4.076055392200811, 2.7563305479657894]
renderView1.CameraFocalPoint = \
  [2.626482081133212, 3.077787159917152, -0.8435425660433065]
renderView1.CameraViewUp = \
  [0.10631851294057819, 0.4261581357409603, 0.898379439406253]
renderView1.CameraParallelScale = \
  1.07460306336
renderView1.CameraParallelProjection = 1


renderView1.Shadows = 0
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.EnableOSPRay = 1
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
surf = LegacyVTKReader(FileNames=ff)

# list of all sources
vs = [surf]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
surf = ForceTime(surf)

# all ForceTime
vft = [surf]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# create a new 'Clip'
if 0:
  clip1 = Clip(Input=surf)
  clip1.ClipType = 'Plane'
  clip1.ClipType.Origin = [0.5, 0.5, 0.5]
  clip1.ClipType.Normal = [0.0, 0.0, 1.0]
  surf = clip1

surfshow = Show(surf, renderView1)
surfshow.Representation = 'Surface'
surfshow.ColorArrayName = ['CELLS', '']
surfshow.Opacity = 0.5

if 0:
  calculator1 = Calculator(Input=surf)
  calculator1.AttributeType = 'Cell Data'
  calculator1.ResultArrayName = 'scl'
  calculator1.Function = 'sin(cl*1234567)'
  calculator1Display = Show(calculator1, renderView1)
  sclLUT = GetColorTransferFunction('scl')
  sclLUT.RGBPoints = [-0.999909626352117, 0.278431372549, 0.278431372549, 0.858823529412, -0.7172826127938886, 0.0, 0.0, 0.360784313725, -0.4366320119178857, 0.0, 1.0, 1.0, -0.15202858567743205, 0.0, 0.501960784314, 0.0, 0.12862201519857086, 1.0, 1.0, 0.0, 0.41124902875679914, 1.0, 0.380392156863, 0.0, 0.6938760423150276, 0.419607843137, 0.0, 0.0, 0.9765030558732559, 0.878431372549, 0.301960784314, 0.301960784314]
  sclLUT.ColorSpace = 'RGB'
  sclLUT.ScalarRangeInitialized = 1.0
  calculator1Display.Representation = 'Surface'
  calculator1Display.ColorArrayName = ['CELLS', 'scl']
  calculator1Display.LookupTable = sclLUT
  calculator1Display.Opacity = 0.5



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
