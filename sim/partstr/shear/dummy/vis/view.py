#!/usr/bin/env pvbatch

from paraview.simple import *

from glob import glob
import sys
import re
import math
import os
import numpy as np

def Error(msg):
    sys.stderr.write(msg + "\n")
    exit(1)

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
    sys.stderr.write('''usage: {:} [vf_*.vtk]
Plots bubbles.
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

# vf input
ff = natsorted(av[1:])
if not len(ff):
    Error("empty file list")
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

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1480, 1000]
renderView1.OrientationAxesVisibility = 0
renderView1.CameraPosition = [-0.390147852074483, 1.9793716269924935, 4.31939983327799]
renderView1.CameraFocalPoint = [1.0, 0.5000071823596954, 0.5]
renderView1.CameraViewUp = [0.11697777844051119, 0.9396926207859085, -0.32139380484326946]
renderView1.CameraParallelScale = 0.7551670624798756
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.0]*3
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5

# ----------------------------------------------------------------
# END RENDER
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
ext = os.path.splitext(ff[0])[1]
if ext == ".vtk":
    vf = LegacyVTKReader(FileNames=ff)
elif ext == ".xmf":
    vf = XDMFReader(FileNames=ff)
    vf.CellArrayStatus = ['vf']
    vf.GridStatus = ['Grid_0']
else:
    Error("unknown extension '{:}'".format(ext))


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

# create a new 'Cell Data to Point Data'
cellpnt = CellDatatoPointData(Input=vf)
cellpnt.CellDataArraytoprocess = ['f']

cont = Contour(Input=cellpnt)
cont.ContourBy = ['POINTS', 'f']
cont.Isosurfaces = [0.5]
cont.PointMergeMethod = 'Uniform Binning'

isovol = IsoVolume(Input=cellpnt)
isovol.InputScalars = ['POINTS', 'f']
isovol.ThresholdRange = [0.5, 1.1]

extr = ExtractSurface(Input=isovol)

generateSurfaceNormals1 = GenerateSurfaceNormals(Input=extr)

clip = Clip(Input=generateSurfaceNormals1)
clip.ClipType = 'Box'
clip.Scalars = ['POINTS', 'f']
clip.Invert = 0

clip.ClipType.Position = [0.0001, 0.0001, 0.0001]
clip.ClipType.Length = [1.999, 0.999, 0.999]

appnd = AppendDatasets(Input=[cont, clip])

appndDisplay = Show(appnd, renderView1)
appndDisplay.Representation = 'Surface'
appndDisplay.ColorArrayName = [None, '']


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

exit(0)
