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
renderView1.ViewSize = [1000, 2000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 2, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [8.535405181306203, -11.84729921683567, 8.958582336091343]
renderView1.CameraFocalPoint = [1.7784384184717228, 6.717314388057494, 1.7679684977243724]
renderView1.CameraViewUp = [-0.11697777844051055, 0.3213938048432693, 0.9396926207859086]
renderView1.CameraParallelScale = 4.596149164737374
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.AmbientSamples = 10
renderView1.SamplesPerPixel = 10
renderView1.Shadows = 1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelFontFile = ''

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
clpnt = CellDatatoPointData(Input=vf)

# create a new 'Contour'
contour1 = Contour(Input=clpnt)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from clpnt
#clpntDisplay = Show(clpnt, renderView1)
#clpntDisplay.Representation = 'Outline'
#clpntDisplay.ColorArrayName = ['POINTS', '']
#clpntDisplay.AmbientColor = [1.0, 1.0, 1.0]

# show data from contour1
contour1Display = Show(contour1, renderView1)
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']


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
    SaveScreenshot(fn, renderView1, ImageResolution=[2000,2000])

exit(0)
