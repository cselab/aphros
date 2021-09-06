#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np

import paratools

parser = argparse.ArgumentParser(
    description="Renders vorticity magnitude with LAMMPS polymers")
parser.add_argument('--lammps', nargs='+', help="List of data files 'a_*.vtk'")
parser.add_argument('--aphros',
                    nargs='+',
                    help="List of data files 'omm_*.xmf'")
parser.add_argument('--force',
                    action="store_true",
                    help="Overwrite existing files")
parser.add_argument('--draft',
                    action="store_true",
                    help="Fewer samples and lower resolution")
parser.add_argument('--res',
                    type=list,
                    default=[1080, 1080],
                    help="Image resolution in pixels")
parser.add_argument('--samples',
                    type=int,
                    default=10,
                    help="Number of samples per pixel")
parser.add_argument('--colormap', type=str, default="blue_yellow_red")
args = parser.parse_args()

sources_ft = []
timearrays = []

files_omm = args.aphros
source_omm = XDMFReader(FileNames=files_omm)
source_omm.CellArrayStatus = ['omm']
source_omm.GridStatus = ['Grid_10220']
(source_omm, ), (timearray, ) = paratools.ApplyForceTime([source_omm])
sources_ft.append(source_omm)
timearrays.append(timearray)

files_polymer = args.lammps
source_polymer = LegacyVTKReader(FileNames=files_polymer)
(source_polymer, ), (timearray, ) = paratools.ApplyForceTime([source_polymer])
sources_ft.append(source_polymer)
timearrays.append(timearray)

viewsize = args.res
if args.draft:
    viewsize = [max(1, s // 2) for s in viewsize]
    args.samples = min(10, args.samples // 5)

renderView1 = CreateView('RenderView')
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.CameraPosition = [
    3.2933362831074464, 2.111426020928463, 1.6757588679166233
]
renderView1.CameraFocalPoint = [
    0.502268266864121, 0.5000021504238248, 0.5027382206171751
]
renderView1.CameraViewUp = [
    -0.296198132726024, -0.17101007166283444, 0.9396926207859084
]
renderView1.CameraParallelScale = 0.709257758997805
renderView1.CameraParallelProjection = 1
renderView1.ViewSize = viewsize
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.Background = [1] * 3
renderView1.SamplesPerPixel = 1 if args.draft else args.samples

tube1 = Tube(registrationName='Tube1', Input=source_polymer)
tube1.Scalars = ['POINTS', '']
tube1.Vectors = ['POINTS', '1']
tube1.Radius = 0.005

glyph1 = Glyph(registrationName='Glyph1',
               Input=source_polymer,
               GlyphType='Sphere')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.GlyphMode = 'All Points'
glyph1.GlyphType.Radius = 0.1
glyph1.ScaleFactor = 0.09989973691408523

tube1Display = Show(tube1, renderView1, 'GeometryRepresentation')
tube1Display.Representation = 'Surface'
tube1Display.AmbientColor = [0.0, 0.0, 0.0]
tube1Display.ColorArrayName = [None, '']
tube1Display.DiffuseColor = [0.0, 0.0, 0.0]
tube1Display.Opacity = 0.5
tube1Display.Ambient = 0.25

glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']

ommDisplay = Show(source_omm, renderView1, 'UniformGridRepresentation')
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [
    3.0, 0.0, 1.0, 1.0, 7.05, 0.0, 0.0, 1.0, 7.5, 0.0, 0.0, 0.501960784314,
    7.95, 1.0, 0.0, 0.0, 12.0, 1.0, 1.0, 0.0
]
ommLUT.ColorSpace = 'RGB'
ommLUT.ScalarRangeInitialized = 1.0
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.Points = [3.0, 0.0, 0.5, 0.0, 12.0, 1.0, 0.5, 0.0]
ommPWF.ScalarRangeInitialized = 1
ommDisplay.Representation = 'Volume'
ommDisplay.ColorArrayName = ['CELLS', 'omm']
ommDisplay.LookupTable = ommLUT
ommDisplay.SetScaleArray = [None, '']
ommDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ommDisplay.OpacityArray = [None, '']
ommDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ommDisplay.ScalarOpacityUnitDistance = 0.1
ommDisplay.ScalarOpacityFunction = ommPWF
ommDisplay.OpacityArrayName = ['CELLS', 'omm']
ommDisplay.Shade = 1

# FIXME: workaround, otherwise `omm` not shown in the first image
paratools.SetTimeStep(1, sources_ft, timearrays)

pattern = "a_{}.png"
steps = paratools.GetSteps(args.aphros)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force,
                        pattern=pattern)
