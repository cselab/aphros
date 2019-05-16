#!/usr/bin/env pvbatch

# state file generated using paraview version 5.6.0


import sys
import os
import glob
import re
import numpy as np


def F(f):
    return os.path.join(dat, f)

def FF(f):
    return sorted(glob.glob(os.path.join(dat, f)))

def Log(s):
    s += "\n"
    o = sys.stderr
    o.write(s)
    o.flush()

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

def GetCenter(o):
    o.UpdatePipeline()
    di = o.GetDataInformation()
    lim = di.DataInformation.GetBounds()
    lim0 = np.array(lim[::2])
    lim1 = np.array(lim[1::2])
    return (lim0 + lim1) * 0.5

def Norm(v):
    return sum(v ** 2) ** 0.5

def A(v):
    return np.array(v)

# rv: render view
# o: object
# bx0, bx1: box
# bxc: center
# cp,cf,cv: reference position, focalpoint and viewup
def ZoomToBox(rv, bx0, bx1, bxc, cp, cf, cv):
    # scale
    sc = np.max(bx1 - bx0) * 1.4
    # normal
    ln = cp - cf
    ln *= sc / Norm(ln)

    k = 1
    rv.CameraPosition = cp * (1 - k) +  bxc * k + ln * 2
    rv.CameraFocalPoint = cf * (1 - k) + bxc * k + ln  * 1
    rv.CameraViewUp = cv

av = sys.argv
if len(av) < 3:
    sys.stderr.write('''usage: {:} dat out
dat: folder with s_*.vtk
out: folder to output png
'''.format(av[0]))
    exit(1)


# data folder
dat = os.path.abspath(av[1])
# output
out = os.path.abspath(av[2])

# parameters
pinter = True
ppart = True

ff = FF("s_*.vtk")
ffb = list(map(os.path.basename, ff))
ffd = list(map(os.path.dirname, ff))
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]

# output pattern (:0 substituted by frame number)
bo = out + "/a_{:}.png"

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000, 500]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.0, 0.2896459652110934, 0.0]
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.StereoType = 0
renderView1.CameraPosition = [0.0, 1.2551414045779636, 5.450255993813145]
renderView1.CameraFocalPoint = [0.0, 0.2852455921685893, -0.05029649402748549]
renderView1.CameraViewUp = [0.0, 0.9848077530122081, -0.17364817766693044]
renderView1.CameraParallelScale = 0.6424932460550369
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableOSPRay = 1
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 10
renderView1.LightScale = 0.75
materialLibrary1 = GetMaterialLibrary()
renderView1.OSPRayMaterialLibrary = materialLibrary1


# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

dat_inter = LegacyVTKReader(FileNames=FF('s_*.vtk'))

# list of all sources
vs = [dat_inter]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
dat_inter = ForceTime(dat_inter)

# all ForceTime
vft = [dat_inter]

# create a new 'Box'
box1 = Box()
box1.XLength = 2.0
box1.YLength = 0.01
box1.ZLength = 2.0
box1.Center = [0.0, -0.005, 0.0]

# create a new 'Transform'
transform1 = Transform(Input=dat_inter)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Scale = [-1.0, 1.0, 1.0]

# create a new 'Append Geometry'
appendGeometry1 = AppendGeometry(Input=[dat_inter, transform1])

# create a new 'Transform'
transform2 = Transform(Input=appendGeometry1)
transform2.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform2.Transform.Scale = [1.0, 1.0, -1.0]

# create a new 'Append Geometry'
appendGeometry2 = AppendGeometry(Input=[transform2, appendGeometry1])

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from box1
box1Display = Show(box1, renderView1)

# trace defaults for the display properties.
box1Display.Representation = 'Surface'
box1Display.AmbientColor = [0.0, 0.0, 0.0]
box1Display.ColorArrayName = [None, '']

appendGeometry2Display = Show(appendGeometry2, renderView1)
appendGeometry2Display.Representation = 'Surface'
appendGeometry2Display.DiffuseColor = np.array([116, 190, 249]) / 255.
appendGeometry2Display.ColorArrayName = [None, '']

SetActiveSource(None)

for i in list(range(len(ss))):
    SetTime(i)
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue
    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1, CompressionLevel=9)
