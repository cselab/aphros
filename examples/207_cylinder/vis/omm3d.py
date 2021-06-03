#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np

import paratools

parser = argparse.ArgumentParser(
    description="Renders geometry and vorticity magnitude.")
parser.add_argument('files', nargs='+', help="List of data files 'omm_*.xmf'")
parser.add_argument('--force',
                    action="store_true",
                    help="Overwrite existing files")
parser.add_argument('--draft', action="store_true", help="less samples")
parser.add_argument('--ebvf',
                    action="store_true",
                    help="Multiply fields by ebvf to hide values in cut cells")
parser.add_argument('--omm_max',
                    type=float,
                    default=150,
                    help="Maximum vorticity magnitude omm, range will be [0,max]")
parser.add_argument('--field_name',
                    type=str,
                    default='omm',
                    help="Name of scalar field for volume rendering")
parser.add_argument('--subset',
                    type=int,
                    nargs=3,
                    default=[1020, 512, 512],
                    help="Subset to extract")
parser.add_argument('--ray', action="store_true", help="raytracing")
parser.add_argument('--res',
                    type=int,
                    default=1920,
                    help="Image width in pixels")
args = parser.parse_args()

source_bcvtk = LegacyVTKReader(
    FileNames=paratools.ReplaceFilename(args.files[:1], "bc.vtk"))

sources_ft = []
timearrays = []

div = 8
viewsize = [1920, 1080]
viewsize = [
    int(s * args.res / viewsize[0] + div - 1) // div * div for s in viewsize
]
if args.field_name != 'omm':
    source_omz = XDMFReader(FileNames=args.files)
    source_omz.CellArrayStatus = [args.field_name]
    source_omz.GridStatus = ['Grid_10220']
    (source_omz, ), (timearray, ) = paratools.ApplyForceTime([source_omz])
    sources_ft.append(source_omz)
    timearrays.append(timearray)
    source_omm = Calculator(Input=source_omz)
    source_omm.AttributeType = 'Cell Data'
    source_omm.ResultArrayName = 'omm'
    source_omm.Function = 'abs({})'.format(args.field_name)
else:
    source_omm = XDMFReader(FileNames=args.files)
    source_omm.CellArrayStatus = ['omm']
    source_omm.GridStatus = ['Grid_10220']
    (source_omm, ), (timearray, ) = paratools.ApplyForceTime([source_omm])
    sources_ft.append(source_omm)
    timearrays.append(timearray)

if args.ebvf:
    files_ebvf = paratools.ReplaceFilename(args.files, "ebvf_{}.xmf")
    source_ebvf = XDMFReader(FileNames=files_ebvf)
    source_ebvf.CellArrayStatus = ['ebvf']
    source_ebvf.GridStatus = ['Grid_10220']
    (source_ebvf, ), (timearray, ) = paratools.ApplyForceTime([source_ebvf])
    sources_ft.append(source_ebvf)
    timearrays.append(timearray)


if args.ebvf:
    source_omm = AppendAttributes(Input=[source_omm, source_ebvf])
    calc_omm = Calculator(Input=source_omm)
    calc_omm.AttributeType = 'Cell Data'
    calc_omm.ResultArrayName = 'omm'
    calc_omm.Function = 'min(1, omm / {:}) * (ebvf^4)'.format(args.omm_max)
else:
    calc_omm = Calculator(Input=source_omm)
    calc_omm.AttributeType = 'Cell Data'
    calc_omm.ResultArrayName = 'omm'
    calc_omm.Function = 'min(1, omm / {:})'.format(args.omm_max)

calc_omm = ExtractSubset(Input=calc_omm)
calc_omm.VOI = [0, args.subset[0], 0, args.subset[1], 0, args.subset[2]]

renderView1 = CreateView('RenderView')
renderView1.ViewSize = viewsize
renderView1.OrientationAxesVisibility = 0
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.CameraPosition = [
    0.6434987006577046, 2.6590788377698193, 4.67979745197099
]
renderView1.CameraFocalPoint = [
    6.963243379949931, -2.2383662155954815, -7.089237203804379
]
renderView1.CameraViewUp = [
    0.20574948679114008, 0.9377672763684191, -0.27974932360550886
]
renderView1.Background = [1] * 3
if args.ray:
    renderView1.EnableRayTracing = 1
    renderView1.BackEnd = 'OSPRay raycaster'
    renderView1.SamplesPerPixel = 1 if args.draft else 10

threshold_bcvtk = Threshold(Input=source_bcvtk)
threshold_bcvtk.Scalars = ['CELLS', 'group']
threshold_bcvtk.ThresholdRange = [3, 3]
bcvtkDisplay = Show(threshold_bcvtk, renderView1)
bcvtkDisplay.Representation = 'Surface'
bcvtkDisplay.ColorArrayName = [None, '']
bcvtkDisplay.AmbientColor = [1] * 3
bcvtkDisplay.DiffuseColor = [1] * 3

ommDisplay = Show(calc_omm, renderView1, 'UniformGridRepresentation')
ommLUT = GetColorTransferFunction('omm')
ommLUT.AutomaticRescaleRangeMode = 'Never'
ommLUT.RGBPoints = [
    0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.143, 0.0, 0.0,
    0.360784313725, 0.285, 0.0, 1.0, 1.0, 0.429, 0.0, 0.501960784314, 0.0,
    0.571, 1.0, 1.0, 0.0, 0.714, 1.0, 0.380392156863, 0.0, 0.857,
    0.419607843137, 0.0, 0.0, 1.0, 0.878431372549, 0.301960784314,
    0.301960784314
]
ommLUT.ColorSpace = 'RGB'
ommLUT.ScalarRangeInitialized = 1.0
ommPWF = GetOpacityTransferFunction('omm')
ommPWF.ScalarRangeInitialized = 1
ommDisplay.Representation = 'Volume'
ommDisplay.ColorArrayName = ['CELLS', 'omm']
ommDisplay.LookupTable = ommLUT
ommDisplay.ScalarOpacityUnitDistance = 0.04
ommDisplay.Shade = 1
ommDisplay.ScalarOpacityFunction = ommPWF
ommDisplay.OpacityTransferFunction.Points = [
    0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0
]

steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)

