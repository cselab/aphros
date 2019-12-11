#!/usr/bin/env pvbatch

from paraview.simple import *

from glob import glob
import sys
import re
import math
import os
import numpy as np

def IsInteractive():
    return "__file__" not in globals() and "vtkpython" in sys.executable

interactive = IsInteractive()

def GetArgs():
    if interactive:
        import PyQt5.QtWidgets
        return PyQt5.QtWidgets.QFileDialog.getOpenFileNames(None, "Open input files")[0]

    av = sys.argv
    if len(av) < 2 or av[1] == '-h':
        sys.stderr.write(
'''usage: {:} [sm_*.vtk]
Plots bubbles.
# Output:
# a_*.png in current folder
'''.format(av[0]))
        exit(1)
    return av[1:]



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

av = GetArgs()

# vf input
ff = natsorted(av)
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
renderView1.ViewSize = [1000,1000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 1.5, 0.5]
renderView1.StereoType = 0
renderView1.UseLight = 1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.CameraPosition =\
[2.0641204031229305, 1.3469838800348422, 3.208017847685781]
renderView1.CameraFocalPoint =\
[0.44809499933849006, 0.48095847625040666, 0.40897974200911913]
renderView1.CameraViewUp =\
[-0.12940952255125976, 0.9659258262890686, -0.22414386804201236]
renderView1.CameraParallelScale = 0.6869798893272092
renderView1.CameraParallelProjection = 1

renderView1.Background = [1.]*3

if interactive:
    layout1 = CreateLayout(name='Layout')
    layout1.AssignView(0, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
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

surf = GenerateSurfaceNormals(Input=surf)

surf = Calculator(Input=surf)
surf.ResultNormals = 1
surf.AttributeType = 'Point Data'
surf.ResultArrayName = 'normals'
surf.Function = 'nn'

clip1 = Clip(Input=surf)
clip1.ClipType = 'Plane'
clip1.Crinkleclip = 1
clip1.ClipType.Origin = [0.0, 0.99, 0.0]
clip1.ClipType.Normal = [0.0, 1.0, 0.0]
surf=clip1

surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = 'Surface'
surfDisplay.AmbientColor = [0.0, 0.0, 0.0]
surfDisplay.Opacity = 1
surfDisplay.ColorArrayName = [None, '']

#####################################################
### END OF STATE FILE
#####################################################

if interactive:
    SetActiveSource(sm_)
else:
    for i in list(range(len(ss))):
        fn = bo.format("{:04d}".format(ss[i]))
        if os.path.isfile(fn):
            Log("skip existing {:}".format(fn))
            continue

        SetTime(i)

        Log("{:}/{:}: {:}".format(i + 1, len(ss), fn))
        SaveScreenshot(fn, renderView1)
