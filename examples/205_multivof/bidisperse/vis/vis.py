#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np

import paratools

parser = argparse.ArgumentParser(
    description="Renders geometry, bubbles and vorticity field.")
parser.add_argument('files', nargs='+', help="list of data files 'sm_*.vtk'")
parser.add_argument('--force',
                    action="store_true",
                    help="overwrite existing files")
parser.add_argument('--omz',
                    action="store_true",
                    help="slice of z-component of vorticity")
parser.add_argument('--omz_max',
                    type=float,
                    default=150,
                    help="maximum vorticity omz, range will be [-max,max]")
parser.add_argument('--bubble_slice',
                    action="store_true",
                    help="draw countours of bubble slices instead of surfaces")
parser.add_argument('--white',
                    action="store_true",
                    help="white background")
parser.add_argument('--preset_omz',
                    action="store_true",
                    help="combines --omz --bubble_slice --white")
parser.add_argument('--res',
                    type=int,
                    default=2048,
                    help="pixels per unit length (channel width)")
args = parser.parse_args()

if args.preset_omz:
    args.omz = True
    args.bubble_slice = True
    args.white = True

source_bcvtk = LegacyVTKReader(
    FileNames=paratools.ReplaceFilename(args.files[:1], "bc.vtk"))
source_sm = LegacyVTKReader(FileNames=args.files)
sources_ft, timearrays = paratools.ApplyForceTime([source_sm])
source_sm, = sources_ft

if args.omz:
    files_omz = paratools.ReplaceFilename(args.files, "omz_{}.xmf")
    source_omz = XDMFReader(FileNames=files_omz)
    source_omz.CellArrayStatus = ['omz']
    source_omz.GridStatus = ['Grid_10220']
    (source_omz,), (timearray,) = paratools.ApplyForceTime([source_omz])
    sources_ft.append(source_omz)
    timearrays.append(timearray)

box = paratools.GetBoundingBox(source_bcvtk)
boxc = (box[0] + box[1]) * 0.5
pixel = 1. / args.res
linewidth = args.res / 1024.
div = 16
boxsize = box[1] - box[0] + pixel * div
viewsize = [int(a * args.res + div - 1) // div * div for a in boxsize[:2]]

color_gray = [0.75] * 3

renderView1 = CreateView('RenderView')
renderView1.ViewSize = viewsize
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.CameraPosition = [boxc[0], boxc[1], 10]
renderView1.CameraFocalPoint = [boxc[0], boxc[1],  0]
renderView1.CameraParallelScale = 0.5 * boxsize[1] * 1.05
renderView1.CameraParallelProjection = 1
renderView1.Background = color_gray

try:
    light1 = CreateLight()
    light1.Coords = 'Ambient'
    renderView1.AdditionalLights = light1
except:
    pass

if args.white:
    renderView1.Background = [1]*3
    renderView1.AdditionalLights = []

source_sm = Calculator(Input=source_sm)
source_sm.ResultNormals = 1
try:
    source_sm.AttributeType = 'Point Data'
except:
    pass
source_sm.ResultArrayName = 'normals'
source_sm.Function = 'nn'

if args.omz:
    celltopoint_omz = CellDatatoPointData(Input=source_omz)
    celltopoint_omz.CellDataArraytoprocess = ['omz']
    slice_omz = Slice(Input=celltopoint_omz)
    slice_omz.SliceType = 'Plane'
    slice_omz.HyperTreeGridSlicer = 'Plane'
    slice_omz.SliceType.Origin = boxc
    slice_omz.SliceType.Normal = [0.0, 0.0, 1.0]
    slice_omzDisplay = Show(slice_omz, renderView1)
    omzLUT = GetColorTransferFunction('omz')
    # geo color scheme (blue/red)
    omzLUT.RGBPoints = [
            -args.omz_max, 0.0, 0.6039215686274509, 0.8705882352941177,
            0, 1.0, 1.0, 1.0,
            args.omz_max, 1.0, 0.12156862745098039, 0.3568627450980392]
    omzLUT.Discretize = 0
    omzLUT.ScalarRangeInitialized = 1.0
    slice_omzDisplay.Representation = 'Surface'
    slice_omzDisplay.ColorArrayName = ['POINTS', 'omz']
    slice_omzDisplay.LookupTable = omzLUT

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
    slice_smDisplay.LineWidth = 2.0 * linewidth
else:
    sm_Display = Show(source_sm, renderView1)
    sm_Display.Representation = 'Surface'
    sm_Display.AmbientColor = [0.0, 0.0, 0.0]
    sm_Display.ColorArrayName = ['POINTS', '']
    sm_Display.DiffuseColor = [0.0, 0.0, 0.0]
    sm_Display.Specular = 1.0
    sm_Display.SpecularColor = color_gray
    sm_Display.SpecularPower = 10.
    sm_Display.Diffuse = 0.0

slice1 = Slice(Input=source_bcvtk)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = boxc
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1Display = Show(slice1, renderView1)
slice1Display.Representation = 'Wireframe'
slice1Display.AmbientColor = [0.0, 0.0, 0.0]
slice1Display.ColorArrayName = ['POINTS', '']
slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
slice1Display.LineWidth = 4. * linewidth

steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)
