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

# Sets time of datasets to step i
def SetTime(i):
    global vft, vt
    for j in range(len(vft)):
        s = vft[j]
        if i < len(vt[j]):
            s.ForcedTime = vt[j][i]
        s.UpdatePipeline()

av = sys.argv
if len(av) < 2 or av[1] == '-h':
    sys.stderr.write('''usage: {:} [vf_*.xmf]
# Output:
# a_*.png in current folder
'''.format(av[0]))
    exit(1)

av = av[1:]

# vf input
ff = sorted(av[0:])
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

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [250, 500]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [2.0, 4.0, 2.0]
renderView1.UseLight = 0
renderView1.StereoType = 0
renderView1.CameraPosition = [2.0, 4.0, 35.71737210362048]
renderView1.CameraFocalPoint = [2.0, 4.0, 2.0]
renderView1.CameraParallelScale = 4.0753768202092715
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.OSPRayMaterialLibrary = materialLibrary1



# ----------------------------------------------------------------
# BEGIN READERS
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
read_vf = XDMFReader(FileNames=ff)
read_vf.CellArrayStatus = ['vf']
read_vf.GridStatus = ['Grid_0']

read_omm = XDMFReader(FileNames=ffomm)
read_omm.CellArrayStatus = ['omm']
read_omm.GridStatus = ['Grid_1']

# list of all sources
vs = [read_vf, read_omm]

# time steps
vt = [np.array(s.TimestepValues) for s in vs]

# replace with ForceTime
read_vf = ForceTime(read_vf)
read_omm = ForceTime(read_omm)

# all ForceTime
vft = [read_vf, read_omm]


# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

#####################################################
### BEGIN OF STATE FILE
#####################################################

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Cell Data to Point Data'
cellDatatoPointData2 = CellDatatoPointData(Input=read_vf)

# create a new 'Contour'
contour1 = Contour(Input=cellDatatoPointData2)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Slice'
slice1 = Slice(Input=contour1)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [2.0, 4.0, 2.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# create a new 'Slice'
slice2 = Slice(Input=read_omm)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [2.0, 4.0, 2.0]
slice2.SliceType.Normal = [0.0, 0.0, 1.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice2
slice2Display = Show(slice2, renderView1)

# get color transfer function/color map for 'omm'
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.EnableOpacityMapping = 1
ommLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 5.0, 0.0, 0.0, 0.0]
ommLUT.ColorSpace = 'RGB'
ommLUT.NanColor = [1.0, 0.0, 0.0]
ommLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.ColorArrayName = ['CELLS', 'omm']
slice2Display.LookupTable = ommLUT
slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice2Display.SelectOrientationVectors = 'None'
slice2Display.ScaleFactor = 0.8
slice2Display.SelectScaleArray = 'omm'
slice2Display.GlyphType = 'Arrow'
slice2Display.GlyphTableIndexArray = 'omm'
slice2Display.GaussianRadius = 0.04
slice2Display.SetScaleArray = [None, '']
slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
slice2Display.OpacityArray = [None, '']
slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
slice2Display.DataAxesGrid = 'GridAxesRepresentation'
slice2Display.SelectionCellLabelFontFile = ''
slice2Display.SelectionPointLabelFontFile = ''
slice2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice2Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice2Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice2Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from slice1
slice1Display = Show(slice1, renderView1)

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = [None, '']
slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
slice1Display.PointSize = 30.0
slice1Display.LineWidth = 3.0
slice1Display.OSPRayScaleArray = 'Normals'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.4
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.02
slice1Display.SetScaleArray = ['POINTS', 'Normals']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'Normals']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.SelectionCellLabelFontFile = ''
slice1Display.SelectionPointLabelFontFile = ''
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [0.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1

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

exit(0)
