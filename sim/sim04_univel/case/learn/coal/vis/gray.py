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
    sys.stderr.write('''usage: {:} [sm_*.vtk]
Plots isosurface.
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

def CheckFlag(name):
    if name in av:
        av.remove(name)
        return True
    return False


draft = CheckFlag('-draft')
orange = CheckFlag('-orange')
green = CheckFlag('-green')

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

# Create a new 'Light'
light1 = CreateLight()
light1.Intensity = 0.8
light1.Position = [0.7, 1.0, 1.8819245943648655]
light1.FocalPoint = [0.500005267560482, 0.4999827891588211, 0.4998413771390915]

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.OrientationAxesVisibility = 0
renderView1.CameraPosition =\
[0.5, 0.5, 2]
renderView1.CameraFocalPoint =\
[0.5, 0.5, 0.5]
renderView1.CameraViewUp =\
[0.0, 1.0, 0.0]
renderView1.CameraParallelScale = 0.3556556576709247
renderView1.CameraParallelProjection = 0
renderView1.Background = [1.0]*3
renderView1.UseLight = 0
renderView1.OSPRayMaterialLibrary = materialLibrary1
ospray = 1
if hasattr(renderView1, 'EnableOSPray'):
    renderView1.EnableOSPRay = ospray
    renderView1.OSPRayRenderer = 'raycaster'
if hasattr(renderView1, 'EnableRayTracing'):
    renderView1.EnableRayTracing = ospray
    renderView1.BackEnd = 'raycaster'
    renderView1.Denoise = 1
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 5 if draft else 30
renderView1.AdditionalLights = light1


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


vf = CellDatatoPointData(Input=vf)
surf = Contour(Input=vf)
surf.ContourBy = ['POINTS', 'vf']
surf.Isosurfaces = [0.5]
#surf = GenerateSurfaceNormals(Input=surf)
surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = 'Surface'
surfDisplay.ColorArrayName = ['POINTS', '']
surfDisplay.DiffuseColor = [1]*3
#c_green = [0.169, 0.627, 0.169]
c_green = [0.636, 0.788, 0.604]
#c_orange = [1., 0.494, 0.055]
c_orange = [0.997, 0.698, 0.493]
if green: surfDisplay.DiffuseColor = c_green
if orange: surfDisplay.DiffuseColor = c_orange


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
    SaveScreenshot(fn, renderView1, TransparentBackground=0)

exit(0)
