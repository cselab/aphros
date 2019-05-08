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
dat: folder with s_*.vtk, partit_*.csv
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
renderView1.ViewSize = [500, 500]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5253897160291672, 0.5222961753606796, 0.5137387216091156]
renderView1.UseLight = 1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.LightScale = 0.85
renderView1.CameraPosition = [0.24839119124607448, -0.10623063113324308, 0.3704770501303224]
renderView1.CameraFocalPoint = [0.7212417088961819, 0.6928833886588502, 0.5883919766652376]
renderView1.CameraViewUp = [0.7866517756373992, -0.5447244857787051, 0.2906100797970575]

cam = F("cam.py")
if os.path.isfile(cam):
    print("Reading custom camera position")
    print(cam)
    ll = open(cam).readlines()
    renderView1.CameraPosition = eval(ll[0])
    renderView1.CameraFocalPoint = eval(ll[1])
    renderView1.CameraViewUp = eval(ll[2])

cam_pos = A(renderView1.CameraPosition)
cam_foc = A(renderView1.CameraFocalPoint)
cam_view = A(renderView1.CameraViewUp)

renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableOSPRay = 1
renderView1.AmbientSamples = 1
renderView1.SamplesPerPixel = 50
renderView1.ProgressivePasses = 1
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
dat_partcon = LegacyVTKReader(FileNames=FF('partit_*.vtk'))
dat_part = CSVReader(FileName=FF('partit_*.csv'))

# list of all sources
vs = [dat_inter, dat_partcon, dat_part]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
dat_inter = ForceTime(dat_inter)
dat_partcon = ForceTime(dat_partcon)
dat_part = ForceTime(dat_part)

# all ForceTime
vft = [dat_inter, dat_partcon, dat_part]

# create a new 'Programmable Filter'
prog_partcon = ProgrammableFilter(Input=dat_partcon)
prog_partcon.Script = """execfile("prog_part.py") """
prog_partcon.RequestInformationScript = ''
prog_partcon.RequestUpdateExtentScript = ''
prog_partcon.PythonPath = ''

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=dat_part)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Programmable Filter'
prog_part = ProgrammableFilter(Input=tableToPoints1)
prog_part.Script = """execfile("prog_part.py")"""
prog_part.RequestInformationScript = ''
prog_part.RequestUpdateExtentScript = ''
prog_part.PythonPath = ''

# create a new 'Clip'
clippart = Clip(Input=prog_part)
clippart.ClipType = 'Scalar'
clippart.Scalars = ['POINTS', 'sp']
clippart.Value = 0
clippart.Invert = 0

# create a new 'Calculator'
calcpart = Calculator(Input=clippart)
calcpart.AttributeType = 'Point Data'
calcpart.ResultArrayName = 'cc'
calcpart.Function = 'sin(1234.5678*c)'

# create a new 'Clip'
clip_partcon = Clip(Input=prog_partcon)
clip_partcon.ClipType = 'Scalar'
clip_partcon.Scalars = ['CELLS', 'sc']
clip_partcon.Value = 0
clip_partcon.Invert = 0

calcpartcon = Calculator(Input=clip_partcon)
calcpartcon.AttributeType = 'Cell Data'
calcpartcon.ResultArrayName = 'cc'
calcpartcon.Function = 'sin(1234.5678*c)'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# get color transfer function/color map for 'cc'
ccLUT = GetColorTransferFunction('cc')
ccLUT.AutomaticRescaleRangeMode = 'Never'
ccLUT.RGBPoints = [-1.0, 0.0, 0.0, 1.0, -0.6679999999999999, 0.0, 0.0, 1.0, -0.6659999999999999, 1.0, 0.0, 1.0, -0.33599999999999997, 1.0, 0.0, 1.0, -0.33399999999999996, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0020000000000000018, 0.0, 1.0, 0.0, 0.3320000000000001, 0.0, 1.0, 0.0, 0.3340000000000001, 1.0, 1.0, 0.0, 0.6639999999999997, 1.0, 1.0, 0.0, 0.6659999999999999, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
ccLUT.ColorSpace = 'HSV'
ccLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'cc'
ccPWF = GetOpacityTransferFunction('cc')
ccPWF.Points = [-1.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
ccPWF.ScalarRangeInitialized = 1

if pinter:
    interDisplay = Show(dat_inter, renderView1)
    interDisplay.Representation = 'Surface With Edges'
    interDisplay.AmbientColor = [0.0, 0.0, 0.0]
    interDisplay.ColorArrayName = ['POINTS', '']
    interDisplay.DiffuseColor = [0.75, 0.75, 0.75]
    interDisplay.LineWidth = 1
    interDisplay.RenderLinesAsTubes = 1
    interDisplay.EdgeColor = [0.0, 0.0, 0.0]

if ppart:
    clippartDisplay = Show(calcpart, renderView1)
    clippartDisplay.Representation = 'Surface'
    clippartDisplay.AmbientColor = [0.0, 0.0, 0.0]
    clippartDisplay.ColorArrayName = ['POINTS', 'cc']
    clippartDisplay.PointSize = 10.0
    clippartDisplay.RenderPointsAsSpheres = 1
    clippartDisplay.LookupTable = ccLUT

    clip_partconDisplay = Show(calcpartcon, renderView1)
    clip_partconDisplay.Representation = 'Surface'
    clip_partconDisplay.AmbientColor = [1.0, 0.0, 0.0]
    clip_partconDisplay.ColorArrayName = ['CELLS', 'cc']
    clip_partconDisplay.LineWidth = 5
    clip_partconDisplay.RenderLinesAsTubes = 1
    clip_partconDisplay.LookupTable = ccLUT

SetActiveSource(None)

#anim = GetAnimationScene()
#anim.NumberOfFrames = len(dat_inter.TimestepValues)
#anim.PlayMode = 'Snap To TimeSteps'
#SaveAnimation(out, FrameWindow=[0,10], CompressionLevel=9)
#SaveAnimation(out, CompressionLevel=9)

bx0, bx1 = GetBox(dat_inter)

for i in list(range(len(ss))):
    SetTime(i)
    bxc = GetCenter(dat_inter)
    ZoomToBox(renderView1, bx0, bx1, bxc, cam_pos, cam_foc, cam_view)
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue
    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1, CompressionLevel=9)
