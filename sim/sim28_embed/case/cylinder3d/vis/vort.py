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
    sys.stderr.write('''usage: {:} [omm_*.xmf]
Plots isosurface.
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(0)

def CheckFlag(name):
    if name in av:
        av.remove(name)
        return True
    return False

def CheckVar(name):
    if name in av:
        i = av.index(name)
        r = av[i + 1]
        del av[i]
        del av[i]
        return r
    return None


cam = 1  # view from side perspective
if CheckFlag('-C1'):
    cam = 1
if CheckFlag('-C2'):
    cam = 2
draft = CheckFlag('-draft')
vortk = CheckVar('-vortk')
surf = CheckFlag('-surf')

# omm input
ff = natsorted(av[1:])
# sm basename
ffb = list(map(os.path.basename, ff))
# sm dirname
ffd = list(map(os.path.dirname, ff))
# steps
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]
# vf input
ffvf = [os.path.join(d, "ebvf_{:04d}.xmf".format(s)) for d,s in zip(ffd,ss)]
# omm input
ffomm = [os.path.join(d, "omm_{:04d}.xmf".format(s)) for d,s in zip(ffd,ss)]
# s input
ffs = [os.path.join(d, "s_{:04d}.vtk".format(s)) for d,s in zip(ffd,ss)]


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

materialLibrary1 = GetMaterialLibrary()


# view
C1 = [
        [8.023651268489624, -10.053402877533147, 5.531727456729329],
        [1.401771596809996, 16.68028345204133, -2.2833261546897403],
        [-0.04793261079752272, 0.2693077225580858, 0.9618606008111105],
    ]

C2 = [
        [12.832352004532979, -0.7993424868513054, 5.205722861235075],
        [-10.728365828708178, 12.312741594120487, -4.416135812637067],
        [-0.29958171881826634, 0.1527692731533269, 0.941760236434995],
    ]


CC = [C1, C2]
C = CC[cam - 1]

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920,1080]
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 0
renderView1.CameraPosition = C[0]
renderView1.CameraFocalPoint = C[1]
renderView1.CameraViewUp = C[2]

renderView1.Background = [1.0]*3
renderView1.UseLight = 1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

vf = XDMFReader(FileNames=ffvf)
vf.CellArrayStatus = ['ebvf']
vf.GridStatus = ['Grid_0']

omm = XDMFReader(FileNames=ffomm)
omm.CellArrayStatus = ['omm']
omm.GridStatus = ['Grid_1']

bcvtk = LegacyVTKReader(FileNames=['../bc.vtk'])

if surf:
    ssurf = LegacyVTKReader(FileNames=ffs)

# list of all sources
vs = [vf, omm]
if surf:
    vs.append(ssurf)

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
vf = ForceTime(vf)
omm = ForceTime(omm)

# all ForceTime
vft = [vf, omm]

if surf:
    ssurf = ForceTime(ssurf)
    vft.append(ssurf)


# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

appendAttributes1 = AppendAttributes(Input=[omm, vf])
ommeb = Calculator(Input=appendAttributes1)
ommeb.AttributeType = 'Cell Data'
ommeb.ResultArrayName = 'ommeb'
ommeb.Function = 'omm*ebvf'
ommebDisplay = Show(ommeb, renderView1)

ommebLUT = GetColorTransferFunction('ommeb')
ommebLUT.AutomaticRescaleRangeMode = 'Never'
ommebLUT.RGBPoints = [4.0, 0.0, 0.0, 0.34902, 4.5, 0.039216, 0.062745, 0.380392, 5.0, 0.062745, 0.117647, 0.411765, 5.5, 0.090196, 0.184314, 0.45098, 6.0, 0.12549, 0.262745, 0.501961, 6.5, 0.160784, 0.337255, 0.541176, 7.0, 0.2, 0.396078, 0.568627, 7.5, 0.239216, 0.454902, 0.6, 8.0, 0.286275, 0.521569, 0.65098, 8.5, 0.337255, 0.592157, 0.701961, 9.0, 0.388235, 0.654902, 0.74902, 9.5, 0.466667, 0.737255, 0.819608, 10.0, 0.572549, 0.819608, 0.878431, 10.5, 0.654902, 0.866667, 0.909804, 11.0, 0.752941, 0.917647, 0.941176, 11.5, 0.823529, 0.956863, 0.968627, 12.0, 0.941176, 0.984314, 0.988235, 12.0, 0.988235, 0.960784, 0.901961, 12.32, 0.988235, 0.945098, 0.85098, 12.640000000000002, 0.980392, 0.898039, 0.784314, 13.0, 0.968627, 0.835294, 0.698039, 13.5, 0.94902, 0.733333, 0.588235, 14.0, 0.929412, 0.65098, 0.509804, 14.5, 0.909804, 0.564706, 0.435294, 15.0, 0.878431, 0.458824, 0.352941, 15.5, 0.839216, 0.388235, 0.286275, 16.0, 0.760784, 0.294118, 0.211765, 16.5, 0.701961, 0.211765, 0.168627, 17.0, 0.65098, 0.156863, 0.129412, 17.5, 0.6, 0.094118, 0.094118, 18.0, 0.54902, 0.066667, 0.098039, 18.5, 0.501961, 0.05098, 0.12549, 19.0, 0.45098, 0.054902, 0.172549, 19.5, 0.4, 0.054902, 0.192157, 20.0, 0.34902, 0.070588, 0.211765]
ommebLUT.ColorSpace = 'Lab'
ommebLUT.NanColor = [0.25, 0.0, 0.0]
ommebLUT.ScalarRangeInitialized = 1.0
ommebPWF = GetOpacityTransferFunction('ommeb')
ommebPWF.Points = [4.0, 0.0, 0.5, 0.0, 20.0, 1.0, 0.5, 0.0]
ommebPWF.ScalarRangeInitialized = 1

def Mul(v, k):
    for i in range(0, len(v), 4):
        v[i] *= k
if vortk is not None:
    vortk = float(vortk)
    print("using vortk={:}".format(vortk))
    Mul(ommebLUT.RGBPoints, vortk)
    Mul(ommebPWF.Points, vortk)

ommebDisplay.Representation = 'Volume'
ommebDisplay.AmbientColor = [0.0, 0.0, 0.0]
ommebDisplay.ColorArrayName = ['CELLS', 'ommeb']
ommebDisplay.LookupTable = ommebLUT
ommebDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ommebDisplay.SelectOrientationVectors = 'None'
ommebDisplay.OpacityArray = [None, '']
ommebDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ommebDisplay.ScalarOpacityUnitDistance = 0.23
ommebDisplay.ScalarOpacityFunction = ommebPWF
ommebDisplay.Shade = 1

threshold1 = Threshold(Input=bcvtk)
threshold1.Scalars = ['CELLS', 'group']
threshold1Display = Show(threshold1, renderView1)
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['POINTS', '']

def rgb(r, g, b):
    m = 255.
    return [r/m, g/m, b/m]

if surf:
    bubblesDisplay = Show(ssurf, renderView1)
    bubblesDisplay.Representation = 'Surface'
    bubblesDisplay.ColorArrayName = ['POINTS', '']
    bubblesDisplay.DiffuseColor = rgb(255, 127, 14)
    bubblesDisplay.Ambient = 0.2

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
