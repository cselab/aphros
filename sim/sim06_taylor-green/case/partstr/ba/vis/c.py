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


ospray = 0

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
renderView1.EnableOSPRay = ospray
renderView1.AmbientSamples = 10
renderView1.SamplesPerPixel = 10
renderView1.Shadows = 1
renderView1.OSPRayMaterialLibrary = materialLibrary1

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
vf = LegacyVTKReader(FileNames=ff)

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
clpntDisplay = Show(clpnt, renderView1)

# trace defaults for the display properties.
clpntDisplay.Representation = 'Outline'
clpntDisplay.ColorArrayName = ['POINTS', '']
clpntDisplay.OSPRayScaleArray = 'f'
clpntDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clpntDisplay.SelectOrientationVectors = 'None'
clpntDisplay.ScaleFactor = 0.30000000000000004
clpntDisplay.SelectScaleArray = 'f'
clpntDisplay.GlyphType = 'Arrow'
clpntDisplay.GlyphTableIndexArray = 'f'
clpntDisplay.GaussianRadius = 0.015
clpntDisplay.SetScaleArray = ['POINTS', 'f']
clpntDisplay.ScaleTransferFunction = 'PiecewiseFunction'
clpntDisplay.OpacityArray = ['POINTS', 'f']
clpntDisplay.OpacityTransferFunction = 'PiecewiseFunction'
clpntDisplay.DataAxesGrid = 'GridAxesRepresentation'
clpntDisplay.SelectionCellLabelFontFile = ''
clpntDisplay.SelectionPointLabelFontFile = ''
clpntDisplay.PolarAxes = 'PolarAxesRepresentation'
clpntDisplay.AmbientColor = [1.0, 1.0, 1.0]

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clpntDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clpntDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clpntDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clpntDisplay.DataAxesGrid.XTitleFontFile = ''
clpntDisplay.DataAxesGrid.YTitleFontFile = ''
clpntDisplay.DataAxesGrid.ZTitleFontFile = ''
clpntDisplay.DataAxesGrid.XLabelFontFile = ''
clpntDisplay.DataAxesGrid.YLabelFontFile = ''
clpntDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clpntDisplay.PolarAxes.PolarAxisTitleFontFile = ''
clpntDisplay.PolarAxes.PolarAxisLabelFontFile = ''
clpntDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
clpntDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']
contour1Display.OSPRayScaleArray = 'Normals'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.2969206996262074
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.GaussianRadius = 0.014846034981310367
contour1Display.SetScaleArray = ['POINTS', 'Normals']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'Normals']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.SelectionCellLabelFontFile = ''
contour1Display.SelectionPointLabelFontFile = ''
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-0.9999940991401672, 1.0, 0.5, 0.0, 0.9999958872795105, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-0.9999940991401672, 1.0, 0.5, 0.0, 0.9999958872795105, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
contour1Display.DataAxesGrid.XTitleFontFile = ''
contour1Display.DataAxesGrid.YTitleFontFile = ''
contour1Display.DataAxesGrid.ZTitleFontFile = ''
contour1Display.DataAxesGrid.XLabelFontFile = ''
contour1Display.DataAxesGrid.YLabelFontFile = ''
contour1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour1Display.PolarAxes.PolarAxisTitleFontFile = ''
contour1Display.PolarAxes.PolarAxisLabelFontFile = ''
contour1Display.PolarAxes.LastRadialAxisTextFontFile = ''
contour1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(clpnt)
# ----------------------------------------------------------------

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
    SaveScreenshot(fn, renderView1, ImageResolution=[2000,2000])

exit(0)
