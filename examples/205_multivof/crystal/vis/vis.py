#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np

import paratools

parser = argparse.ArgumentParser(
    description="Renders geometry, concentration and particles.")
parser.add_argument('files', nargs='+', help="list of data files 'sm_*.vtk'")
parser.add_argument('--bubble_slice',
                    action="store_true",
                    help="draw countours of bubble slices instead of surfaces")
parser.add_argument('--force',
                    action="store_true",
                    help="overwrite existing files")
args = parser.parse_args()


source_bcvtk = LegacyVTKReader(
    FileNames=paratools.ReplaceFilename(args.files[:1], "bc.vtk"))
source_sm = LegacyVTKReader(FileNames=args.files)
sources_ft, timearrays = paratools.ApplyForceTime([source_sm])
source_sm, = sources_ft

box = paratools.GetBoundingBox(source_bcvtk)
boxc = (box[0] + box[1]) * 0.5
res = 300
pixel = 1. / res
div = 16
boxsize = box[1] - box[0] + pixel * div
viewsize = [int(a * res + div - 1) // div * div for a in boxsize[:2]]


color_gray = [0.8] * 3

renderView1 = CreateView('RenderView')
renderView1.ViewSize = viewsize
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.CameraPosition = [boxc[0], boxc[1], 10]
renderView1.CameraFocalPoint = [boxc[0], boxc[1],  0]
renderView1.CameraParallelScale = boxsize[1] * 0.5 * 1.02
renderView1.CameraParallelProjection = 1
renderView1.Background = color_gray

try:
    light1 = CreateLight()
    light1.Coords = 'Ambient'
    renderView1.AdditionalLights = light1
except:
    pass


slice1 = Slice(Input=source_bcvtk)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = boxc
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from source_sm
source_sm = Calculator(Input=source_sm)
source_sm.ResultNormals = 1
try:
    source_sm.AttributeType = 'Point Data'
except:
    pass
source_sm.ResultArrayName = 'normals'
source_sm.Function = 'nn'


if args.bubble_slice:
    slice_sm = Slice(Input=source_sm)
    slice_sm.SliceType = 'Plane'
    slice_sm.SliceOffsetValues = [0.0]
    slice_sm.SliceType.Origin = boxc
    slice_sm.SliceType.Normal = [0.0, 0.0, 1.0]
    slice_smDisplay = Show(slice_sm, renderView1)
    slice_smDisplay.Representation = 'Surface'
    slice_smDisplay.AmbientColor = [0.0, 0.0, 0.0]
    slice_smDisplay.ColorArrayName = ['POINTS', '']
    slice_smDisplay.DiffuseColor = [0.0, 0.0, 0.0]
    slice_smDisplay.LineWidth = 4.0
else:
    smDisplay = Show(source_sm, renderView1)
    smDisplay.Representation = 'Surface'
    smDisplay.AmbientColor = [0.0, 0.0, 0.0]
    smDisplay.ColorArrayName = ['POINTS', '']
    smDisplay.DiffuseColor = [0.0, 0.0, 0.0]
    smDisplay.Specular = 1.0
    smDisplay.SpecularColor = color_gray
    smDisplay.SpecularPower = 20
    smDisplay.Diffuse = 0.0

slice1Display = Show(slice1, renderView1)
slice1Display.Representation = 'Wireframe'
slice1Display.AmbientColor = [0.0, 0.0, 0.0]
slice1Display.ColorArrayName = ['POINTS', '']
slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
slice1Display.LineWidth = 0.01 / pixel

steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)
