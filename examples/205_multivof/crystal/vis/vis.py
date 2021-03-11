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

light1 = CreateLight()
light1.Coords = 'Ambient'

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
renderView1.AdditionalLights = light1


slice1 = Slice(Input=source_bcvtk)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = boxc
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1.HyperTreeGridSlicer.Origin = [1.5499999523162842, 0.532812487334013, 0.04843749850988388]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from source_sm
source_sm = Calculator(Input=source_sm)
source_sm.ResultNormals = 1
source_sm.AttributeType = 'Point Data'
source_sm.ResultArrayName = 'normals'
source_sm.Function = 'nn'

smDisplay = Show(source_sm, renderView1, 'GeometryRepresentation')
smDisplay.Representation = 'Surface'
smDisplay.AmbientColor = [0.0, 0.0, 0.0]
smDisplay.ColorArrayName = ['POINTS', '']
smDisplay.DiffuseColor = [0.0, 0.0, 0.0]
smDisplay.Specular = 1.0
smDisplay.SpecularColor = color_gray
smDisplay.SpecularPower = 20
smDisplay.Diffuse = 0.0

slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
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
