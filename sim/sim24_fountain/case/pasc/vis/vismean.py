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
    sys.stderr.write('''usage: {:} [sm_*.vtk]
Plots surface of volume fraction averaged in z-direction
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

def CheckFlag(name):
    if name in av:
        av.remove(name)
        return True
    return False

# sm input
ff = natsorted(av[1:])
# sm basename
ffb = list(map(os.path.basename, ff))
# sm dirname
ffd = list(map(os.path.dirname, ff))
# steps
ss = [int(re.findall("_([0-9]*)", fb)[0]) for fb in ffb]
print("ss=", ss)

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

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000,500]
renderView1.UseLight = 0
renderView1.CameraPosition = [1.0, 0.5, 3.575289729747054]
renderView1.CameraFocalPoint = [1.0, 0.5, 0.005208333333333333]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.5
renderView1.CameraParallelProjection = 1
renderView1.OrientationAxesVisibility = 0
renderView1.Background = [1.]*3


# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

vf = XDMFReader(FileNames=ff)
vf.CellArrayStatus = ['data']
vf.GridStatus = ['Grid_0']

# list of all sources
vs = [vf]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
vf = ForceTime(vf)

# all ForceTime
vft = [vf]

# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

celltopnt = CellDatatoPointData(Input=vf)
celltopnt.CellDataArraytoprocess = ['vfz']
celltopntDisplay = Show(celltopnt, renderView1)
dataLUT = GetColorTransferFunction('vfz')
dataLUT.RGBPoints = [0.0, 0.031373, 0.188235, 0.419608, 0.062745, 0.031373, 0.253195, 0.516063, 0.12549, 0.031757, 0.318139, 0.612149, 0.1882355, 0.080969, 0.38113, 0.661361, 0.2509805, 0.130427, 0.444152, 0.710327, 0.3137255, 0.195386, 0.509112, 0.743791, 0.3764705, 0.260715, 0.573841, 0.777209, 0.4392155, 0.341423, 0.628958, 0.808704, 0.501960785, 0.422745, 0.684075, 0.839892, 0.564706, 0.523137, 0.739193, 0.861546, 0.627451, 0.622684, 0.793464, 0.883429, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]
dataLUT.ColorSpace = 'Lab'
dataLUT.ScalarRangeInitialized = 1.0
dataPWF = GetOpacityTransferFunction('vfz')
dataPWF.ScalarRangeInitialized = 1
celltopntDisplay.Representation = 'Surface'
celltopntDisplay.ColorArrayName = ['POINTS', 'vfz']
celltopntDisplay.LookupTable = dataLUT


#####################################################
### END OF STATE FILE
#####################################################

SetTime(1)

for i in list(range(len(ss))):
    fn = bo.format("{:04d}".format(ss[i]))
    if os.path.isfile(fn):
        Log("skip existing {:}".format(fn))
        continue

    SetTime(i)

    Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
    SaveScreenshot(fn, renderView1)
