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
renderView1.ViewSize = [1800, 2000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [3.1415927410125732, 3.1415927410125732, 3.1415927410125732]
renderView1.StereoType = 0
renderView1.CameraPosition = [6.43080086005585, -15.166333096764594, 9.914876482661194]
renderView1.CameraFocalPoint = [3.0001989267510476, 4.2895772779323025, 2.724262644294233]
renderView1.CameraViewUp = [-0.05939117461388465, 0.33682408883346476, 0.9396926207859086]
renderView1.CameraParallelScale = 4.2537082476917885
renderView1.CameraParallelProjection = 1
renderView1.Background = [1., 1., 1.]
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
cellDatatoPointData1 = CellDatatoPointData(Input=vf)

# create a new 'Contour'
contour1 = Contour(Input=cellDatatoPointData1)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Calculator'
calculator1 = Calculator(Input=vf)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'vel'
calculator1.Function = 'iHat*u.x+jHat*u.y+kHat*u.z'

# create a new 'Gradient Of Unstructured DataSet'
grad = GradientOfUnstructuredDataSet(Input=calculator1)
grad.ScalarArray = ['CELLS', 'vel']
grad.ComputeGradient = 0
grad.ComputeVorticity = 1
grad.VorticityArrayName = 'om'


# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from vf
vfDisplay = Show(vf, renderView1)

# trace defaults for the display properties.
vfDisplay.Representation = 'Outline'
vfDisplay.ColorArrayName = ['CELLS', '']
vfDisplay.LineWidth = 3.0
vfDisplay.AmbientColor = [0.0, 0.0, 0.0]


# show data from grad
gradDisplay = Show(grad, renderView1)

# get color transfer function/color map for 'omm'
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [1.0, 0.231373, 0.298039, 0.752941, 3.0000000000000004, 0.865003, 0.865003, 0.865003, 5.0, 0.705882, 0.0156863, 0.14902]
ommLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'omm'
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [1.0, 0.0, 0.5, 0.0, 3.1808218955993652, 0.07236842066049576, 0.5, 0.0, 5.0, 0.2565789520740509, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
gradDisplay.Representation = 'Volume'
gradDisplay.ColorArrayName = ['CELLS', 'om']
gradDisplay.LookupTable = ommLUT
gradDisplay.LineWidth = 3.0
gradDisplay.SelectOrientationVectors = 'None'
gradDisplay.ScaleFactor = 0.628318530717957
gradDisplay.SelectScaleArray = 'None'
gradDisplay.GlyphType = 'Arrow'
gradDisplay.GlyphTableIndexArray = 'None'
gradDisplay.GaussianRadius = 0.03141592653589785
gradDisplay.SetScaleArray = ['POINTS', 'omm']
gradDisplay.ScaleTransferFunction = 'PiecewiseFunction'
gradDisplay.OpacityArray = ['POINTS', 'omm']
gradDisplay.OpacityTransferFunction = 'PiecewiseFunction'
gradDisplay.DataAxesGrid = 'GridAxesRepresentation'
gradDisplay.SelectionCellLabelFontFile = ''
gradDisplay.SelectionPointLabelFontFile = ''
gradDisplay.PolarAxes = 'PolarAxesRepresentation'
gradDisplay.ScalarOpacityUnitDistance = 0.04251092259923938
gradDisplay.ScalarOpacityFunction = ommPWF
gradDisplay.Slice = 128

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
gradDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
gradDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
gradDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
gradDisplay.DataAxesGrid.XTitleFontFile = ''
gradDisplay.DataAxesGrid.YTitleFontFile = ''
gradDisplay.DataAxesGrid.ZTitleFontFile = ''
gradDisplay.DataAxesGrid.XLabelFontFile = ''
gradDisplay.DataAxesGrid.YLabelFontFile = ''
gradDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
gradDisplay.PolarAxes.PolarAxisTitleFontFile = ''
gradDisplay.PolarAxes.PolarAxisLabelFontFile = ''
gradDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
gradDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']
contour1Display.OSPRayScaleArray = 'Normals'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.6283185482025146
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.GaussianRadius = 0.031415927410125735
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
contour1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

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
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
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
    SaveScreenshot(fn, renderView1)

exit(0)
