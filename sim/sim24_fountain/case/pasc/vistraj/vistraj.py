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
renderView1.ViewSize = [1500, 1000]
renderView1.OrientationAxesVisibility = 0
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.CameraPosition = [1.2927692121339105, 0.31638904427303083, 4.760182365073278]
renderView1.CameraFocalPoint = [1.2927692121339105, 0.31638904427303083, 0.49581929395400004]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.35626158103213174
renderView1.CameraParallelProjection = 1
renderView1.Background = [1.]*3
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.AmbientSamples = 1
renderView1.OSPRayMaterialLibrary = materialLibrary1

# ----------------------------------------------------------------
# END RENDER
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

smcl1 = LegacyVTKReader(FileNames=[])
smcl2 = LegacyVTKReader(FileNames=[])
smcl3 = LegacyVTKReader(FileNames=[])

group_sm = GroupDatasets(Input=[smcl1, smcl2, smcl3])

extracttimeAB = ExtractTimeSteps(Input=group_sm)
extracttimeAB.TimeStepIndices = [110, 240]
grouptime = GroupTimeSteps(Input=extracttimeAB)

traj1 = LegacyVTKReader(FileNames=['trajr_7377.vtk'])
traj2 = LegacyVTKReader(FileNames=['trajr_8706.vtk'])
traj3 = LegacyVTKReader(FileNames=['trajr_5228.vtk'])

group_traj = GroupDatasets(Input=[traj3, traj1, traj2])

extract_sm1 = ExtractBlock(Input=grouptime)
extract_sm1.BlockIndices = [2, 6]
extract_sm2 = ExtractBlock(Input=grouptime)
extract_sm2.BlockIndices = [3, 7]
extract_sm3 = ExtractBlock(Input=grouptime)
extract_sm3.BlockIndices = [8, 4]

# create a new 'Legacy VTK Reader'
sm = LegacyVTKReader(FileNames=[])


# ----------------------------------------------------------------
# END READERS
# ----------------------------------------------------------------

# create a new 'Clip'
clip_wall = Clip(Input=sm)
clip_wall.ClipType = 'Scalar'
clip_wall.Scalars = ['CELLS', 'cl']
clip_wall.Value = -0.1
clip_wall.Invert = 0

# create a new 'Calculator'
calcNormals = Calculator(Input=clip_wall)
calcNormals.ResultNormals = 1
calcNormals.ResultArrayName = 'normals'
calcNormals.Function = 'nn'

# create a new 'Clip'
clip_box = Clip(Input=calcNormals)
clip_box.ClipType = 'Box'
clip_box.Scalars = ['CELLS', 'c']
clip_box.Value = 191663739.0

# init the 'Box' selected for 'ClipType'
clip_box.ClipType.Position = [0.0, 0.0, 0.34]
clip_box.ClipType.Length = [2.0, 0.8851796984672546, 0.25]

# create a new 'Force Time'
forcetimeA = ForceTime(Input=clip_box)
forcetimeA.ForcedTime = 110.0

# create a new 'Clip'
clip_cylA = Clip(Input=forcetimeA)
clip_cylA.ClipType = 'Cylinder'
clip_cylA.Scalars = ['CELLS', 'c']
clip_cylA.Value = 191663739.0

# init the 'Cylinder' selected for 'ClipType'
clip_cylA.ClipType.Center = [1.55, 0.2, 0.5]
clip_cylA.ClipType.Axis = [0.0, 0.0, 1.0]
clip_cylA.ClipType.Radius = 0.25

# create a new 'Force Time'
forcetimeB = ForceTime(Input=clip_box)
forcetimeB.ForcedTime = 240.0

# create a new 'Clip'
clip_cylB = Clip(Input=forcetimeB)
clip_cylB.ClipType = 'Cylinder'
clip_cylB.Scalars = ['CELLS', 'c']
clip_cylB.Value = 191663739.0

# init the 'Cylinder' selected for 'ClipType'
clip_cylB.ClipType.Center = [1.05, 0.5, 0.5]
clip_cylB.ClipType.Axis = [0.0, 0.0, 1.0]
clip_cylB.ClipType.Radius = 0.25

# create a new 'Clip'
traj_clipn1 = Clip(Input=group_traj)
traj_clipn1.ClipType = 'Scalar'
traj_clipn1.Scalars = ['POINTS', 'n']
traj_clipn1.Value = 110.0
traj_clipn1.Invert = 0

# create a new 'Clip'
traj_clipn2 = Clip(Input=traj_clipn1)
traj_clipn2.ClipType = 'Scalar'
traj_clipn2.Scalars = ['POINTS', 'n']
traj_clipn2.Value = 240.0

# create a new 'Extract Block'
traj_block2 = ExtractBlock(Input=traj_clipn2)
traj_block2.BlockIndices = [2]

# create a new 'Extract Block'
traj_block1 = ExtractBlock(Input=traj_clipn2)
traj_block1.BlockIndices = [1]

# create a new 'Extract Block'
traj_block3 = ExtractBlock(Input=traj_clipn2)
traj_block3.BlockIndices = [3]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from traj_block1
traj_block1Display = Show(traj_block1, renderView1)

# trace defaults for the display properties.
traj_block1Display.Representation = 'Wireframe'
traj_block1Display.AmbientColor = [0.11764705882352941, 0.4627450980392157, 0.6980392156862745]
traj_block1Display.ColorArrayName = ['POINTS', '']
traj_block1Display.DiffuseColor = [0.11764705882352941, 0.4627450980392157, 0.6980392156862745]
traj_block1Display.LineWidth = 5.0
traj_block1Display.RenderLinesAsTubes = 1
traj_block1Display.Ambient = 0.4
traj_block1Display.OSPRayScaleArray = 'n'
traj_block1Display.OSPRayScaleFunction = 'PiecewiseFunction'
traj_block1Display.SelectOrientationVectors = 'None'
traj_block1Display.ScaleFactor = 0.06610955263
traj_block1Display.SelectScaleArray = 'n'
traj_block1Display.GlyphType = 'Arrow'
traj_block1Display.GlyphTableIndexArray = 'n'
traj_block1Display.GaussianRadius = 0.0033054776315
traj_block1Display.SetScaleArray = ['POINTS', 'n']
traj_block1Display.ScaleTransferFunction = 'PiecewiseFunction'
traj_block1Display.OpacityArray = ['POINTS', 'n']
traj_block1Display.OpacityTransferFunction = 'PiecewiseFunction'
traj_block1Display.DataAxesGrid = 'GridAxesRepresentation'
traj_block1Display.PolarAxes = 'PolarAxesRepresentation'
traj_block1Display.ScalarOpacityUnitDistance = 0.15174790308952285
traj_block1Display.ExtractedBlockIndex = 1

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
traj_block1Display.ScaleTransferFunction.Points = [74.0, 0.0, 0.5, 0.0, 240.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
traj_block1Display.OpacityTransferFunction.Points = [74.0, 0.0, 0.5, 0.0, 240.0, 1.0, 0.5, 0.0]

# show data from traj_block2
traj_block2Display = Show(traj_block2, renderView1)

# trace defaults for the display properties.
traj_block2Display.Representation = 'Wireframe'
traj_block2Display.AmbientColor = [1.0, 0.49411764705882355, 0.054901960784313725]
traj_block2Display.ColorArrayName = ['POINTS', '']
traj_block2Display.DiffuseColor = [1.0, 0.49411764705882355, 0.054901960784313725]
traj_block2Display.LineWidth = 5.0
traj_block2Display.RenderLinesAsTubes = 1
traj_block2Display.Ambient = 0.4
traj_block2Display.OSPRayScaleArray = 'n'
traj_block2Display.OSPRayScaleFunction = 'PiecewiseFunction'
traj_block2Display.SelectOrientationVectors = 'None'
traj_block2Display.ScaleFactor = 0.07751745119320001
traj_block2Display.SelectScaleArray = 'n'
traj_block2Display.GlyphType = 'Arrow'
traj_block2Display.GlyphTableIndexArray = 'n'
traj_block2Display.GaussianRadius = 0.0038758725596600005
traj_block2Display.SetScaleArray = ['POINTS', 'n']
traj_block2Display.ScaleTransferFunction = 'PiecewiseFunction'
traj_block2Display.OpacityArray = ['POINTS', 'n']
traj_block2Display.OpacityTransferFunction = 'PiecewiseFunction'
traj_block2Display.DataAxesGrid = 'GridAxesRepresentation'
traj_block2Display.PolarAxes = 'PolarAxesRepresentation'
traj_block2Display.ScalarOpacityUnitDistance = 0.16944716406487006
traj_block2Display.ExtractedBlockIndex = 1

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
traj_block2Display.ScaleTransferFunction.Points = [86.0, 0.0, 0.5, 0.0, 240.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
traj_block2Display.OpacityTransferFunction.Points = [86.0, 0.0, 0.5, 0.0, 240.0, 1.0, 0.5, 0.0]


# show data from traj_block3
traj_block3Display = Show(traj_block3, renderView1)

# trace defaults for the display properties.
traj_block3Display.Representation = 'Wireframe'
traj_block3Display.AmbientColor = [0.16862745098039217, 0.6274509803921569, 0.16862745098039217]
traj_block3Display.ColorArrayName = ['POINTS', '']
traj_block3Display.DiffuseColor = [0.16862745098039217, 0.6274509803921569, 0.16862745098039217]
traj_block3Display.LineWidth = 5.0
traj_block3Display.RenderLinesAsTubes = 1
traj_block3Display.Ambient = 0.4
traj_block3Display.OSPRayScaleArray = 'n'
traj_block3Display.OSPRayScaleFunction = 'PiecewiseFunction'
traj_block3Display.SelectOrientationVectors = 'None'
traj_block3Display.ScaleFactor = 0.0720227124525
traj_block3Display.SelectScaleArray = 'n'
traj_block3Display.GlyphType = 'Arrow'
traj_block3Display.GlyphTableIndexArray = 'n'
traj_block3Display.GaussianRadius = 0.003601135622625
traj_block3Display.SetScaleArray = ['POINTS', 'n']
traj_block3Display.ScaleTransferFunction = 'PiecewiseFunction'
traj_block3Display.OpacityArray = ['POINTS', 'n']
traj_block3Display.OpacityTransferFunction = 'PiecewiseFunction'
traj_block3Display.DataAxesGrid = 'GridAxesRepresentation'
traj_block3Display.PolarAxes = 'PolarAxesRepresentation'
traj_block3Display.ScalarOpacityUnitDistance = 0.16009256536547248
traj_block3Display.ExtractedBlockIndex = 1

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
traj_block3Display.ScaleTransferFunction.Points = [91.0, 0.0, 0.5, 0.0, 240.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
traj_block3Display.OpacityTransferFunction.Points = [91.0, 0.0, 0.5, 0.0, 240.0, 1.0, 0.5, 0.0]

# show data from extract_sm1
extract_sm1Display = Show(extract_sm1, renderView1)

# trace defaults for the display properties.
extract_sm1Display.Representation = 'Surface'
extract_sm1Display.AmbientColor = [1.0, 0.49411764705882355, 0.054901960784313725]
extract_sm1Display.ColorArrayName = ['POINTS', '']
extract_sm1Display.DiffuseColor = [1.0, 0.49411764705882355, 0.054901960784313725]
extract_sm1Display.PointSize = 30.0
extract_sm1Display.LineWidth = 3.0
extract_sm1Display.RenderPointsAsSpheres = 1
extract_sm1Display.Ambient = 0.4
extract_sm1Display.OSPRayScaleArray = 'nn'
extract_sm1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extract_sm1Display.SelectOrientationVectors = 'nn'
extract_sm1Display.ScaleFactor = 0.03977283239364624
extract_sm1Display.SelectScaleArray = 'c'
extract_sm1Display.GlyphType = 'Arrow'
extract_sm1Display.GlyphTableIndexArray = 'c'
extract_sm1Display.GaussianRadius = 0.001988641619682312
extract_sm1Display.SetScaleArray = ['POINTS', 'nn']
extract_sm1Display.ScaleTransferFunction = 'PiecewiseFunction'
extract_sm1Display.OpacityArray = ['POINTS', 'nn']
extract_sm1Display.OpacityTransferFunction = 'PiecewiseFunction'
extract_sm1Display.DataAxesGrid = 'GridAxesRepresentation'
extract_sm1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extract_sm1Display.ScaleTransferFunction.Points = [-0.9344435930252075, 0.0, 0.5, 0.0, 0.8319071531295776, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extract_sm1Display.OpacityTransferFunction.Points = [-0.9344435930252075, 0.0, 0.5, 0.0, 0.8319071531295776, 1.0, 0.5, 0.0]

# show data from extract_sm2
extract_sm2Display = Show(extract_sm2, renderView1)

# trace defaults for the display properties.
extract_sm2Display.Representation = 'Surface'
extract_sm2Display.AmbientColor = [0.16862745098039217, 0.6274509803921569, 0.16862745098039217]
extract_sm2Display.ColorArrayName = ['POINTS', '']
extract_sm2Display.DiffuseColor = [0.16862745098039217, 0.6274509803921569, 0.16862745098039217]
extract_sm2Display.PointSize = 30.0
extract_sm2Display.LineWidth = 3.0
extract_sm2Display.RenderPointsAsSpheres = 1
extract_sm2Display.Ambient = 0.4
extract_sm2Display.OSPRayScaleArray = 'nn'
extract_sm2Display.OSPRayScaleFunction = 'PiecewiseFunction'
extract_sm2Display.SelectOrientationVectors = 'nn'
extract_sm2Display.ScaleFactor = 0.06540139317512512
extract_sm2Display.SelectScaleArray = 'c'
extract_sm2Display.GlyphType = 'Arrow'
extract_sm2Display.GlyphTableIndexArray = 'c'
extract_sm2Display.GaussianRadius = 0.003270069658756256
extract_sm2Display.SetScaleArray = ['POINTS', 'nn']
extract_sm2Display.ScaleTransferFunction = 'PiecewiseFunction'
extract_sm2Display.OpacityArray = ['POINTS', 'nn']
extract_sm2Display.OpacityTransferFunction = 'PiecewiseFunction'
extract_sm2Display.DataAxesGrid = 'GridAxesRepresentation'
extract_sm2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extract_sm2Display.ScaleTransferFunction.Points = [-0.947035014629364, 0.0, 0.5, 0.0, 0.940422773361206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extract_sm2Display.OpacityTransferFunction.Points = [-0.947035014629364, 0.0, 0.5, 0.0, 0.940422773361206, 1.0, 0.5, 0.0]

# show data from extract_sm3
extract_sm3Display = Show(extract_sm3, renderView1)

# trace defaults for the display properties.
extract_sm3Display.Representation = 'Surface'
extract_sm3Display.AmbientColor = [0.11764705882352941, 0.4627450980392157, 0.6980392156862745]
extract_sm3Display.ColorArrayName = ['POINTS', '']
extract_sm3Display.DiffuseColor = [0.11764705882352941, 0.4627450980392157, 0.6980392156862745]
extract_sm3Display.PointSize = 30.0
extract_sm3Display.LineWidth = 3.0
extract_sm3Display.RenderPointsAsSpheres = 1
extract_sm3Display.Ambient = 0.4
extract_sm3Display.OSPRayScaleArray = 'nn'
extract_sm3Display.OSPRayScaleFunction = 'PiecewiseFunction'
extract_sm3Display.SelectOrientationVectors = 'nn'
extract_sm3Display.ScaleFactor = 0.05287367105484009
extract_sm3Display.SelectScaleArray = 'c'
extract_sm3Display.GlyphType = 'Arrow'
extract_sm3Display.GlyphTableIndexArray = 'c'
extract_sm3Display.GaussianRadius = 0.0026436835527420045
extract_sm3Display.SetScaleArray = ['POINTS', 'nn']
extract_sm3Display.ScaleTransferFunction = 'PiecewiseFunction'
extract_sm3Display.OpacityArray = ['POINTS', 'nn']
extract_sm3Display.OpacityTransferFunction = 'PiecewiseFunction'
extract_sm3Display.DataAxesGrid = 'GridAxesRepresentation'
extract_sm3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extract_sm3Display.ScaleTransferFunction.Points = [-0.9215561747550964, 0.0, 0.5, 0.0, 0.7954306602478027, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extract_sm3Display.OpacityTransferFunction.Points = [-0.9215561747550964, 0.0, 0.5, 0.0, 0.7954306602478027, 1.0, 0.5, 0.0]

# show data from clip_cylB
clip_cylBDisplay = Show(clip_cylB, renderView1)

# trace defaults for the display properties.
clip_cylBDisplay.Representation = 'Surface'
clip_cylBDisplay.AmbientColor = [0.0, 0.0, 0.0]
clip_cylBDisplay.ColorArrayName = ['POINTS', '']
clip_cylBDisplay.Opacity = 0.15
clip_cylBDisplay.OSPRayScaleArray = 'nn'
clip_cylBDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clip_cylBDisplay.SelectOrientationVectors = 'nn'
clip_cylBDisplay.ScaleFactor = 0.1688839316368103
clip_cylBDisplay.SelectScaleArray = 'c'
clip_cylBDisplay.GlyphType = 'Arrow'
clip_cylBDisplay.GlyphTableIndexArray = 'c'
clip_cylBDisplay.GaussianRadius = 0.008444196581840516
clip_cylBDisplay.SetScaleArray = ['POINTS', 'nn']
clip_cylBDisplay.ScaleTransferFunction = 'PiecewiseFunction'
clip_cylBDisplay.OpacityArray = ['POINTS', 'nn']
clip_cylBDisplay.OpacityTransferFunction = 'PiecewiseFunction'
clip_cylBDisplay.DataAxesGrid = 'GridAxesRepresentation'
clip_cylBDisplay.PolarAxes = 'PolarAxesRepresentation'
clip_cylBDisplay.ScalarOpacityUnitDistance = 0.01666272593515284

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip_cylBDisplay.ScaleTransferFunction.Points = [-0.9742168188095093, 0.0, 0.5, 0.0, 0.9698622822761536, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip_cylBDisplay.OpacityTransferFunction.Points = [-0.9742168188095093, 0.0, 0.5, 0.0, 0.9698622822761536, 1.0, 0.5, 0.0]

# show data from clip_cylA
clip_cylADisplay = Show(clip_cylA, renderView1)

# trace defaults for the display properties.
clip_cylADisplay.Representation = 'Surface'
clip_cylADisplay.AmbientColor = [0.0, 0.0, 0.0]
clip_cylADisplay.ColorArrayName = ['POINTS', '']
clip_cylADisplay.Opacity = 0.2
clip_cylADisplay.OSPRayScaleArray = 'nn'
clip_cylADisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clip_cylADisplay.SelectOrientationVectors = 'nn'
clip_cylADisplay.ScaleFactor = 0.0992315411567688
clip_cylADisplay.SelectScaleArray = 'c'
clip_cylADisplay.GlyphType = 'Arrow'
clip_cylADisplay.GlyphTableIndexArray = 'c'
clip_cylADisplay.GaussianRadius = 0.00496157705783844
clip_cylADisplay.SetScaleArray = ['POINTS', 'nn']
clip_cylADisplay.ScaleTransferFunction = 'PiecewiseFunction'
clip_cylADisplay.OpacityArray = ['POINTS', 'nn']
clip_cylADisplay.OpacityTransferFunction = 'PiecewiseFunction'
clip_cylADisplay.DataAxesGrid = 'GridAxesRepresentation'
clip_cylADisplay.PolarAxes = 'PolarAxesRepresentation'
clip_cylADisplay.ScalarOpacityUnitDistance = 0.01298734893503781

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip_cylADisplay.ScaleTransferFunction.Points = [-0.9775257110595703, 0.0, 0.5, 0.0, 0.9647583365440369, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip_cylADisplay.OpacityTransferFunction.Points = [-0.9775257110595703, 0.0, 0.5, 0.0, 0.9647583365440369, 1.0, 0.5, 0.0]


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
