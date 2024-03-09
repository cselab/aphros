#!/usr/bin/env pvbatch

# state file generated using paraview version 5.9.0

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
import argparse
import paratools
import numpy as np
import re
import os

parser = argparse.ArgumentParser(
    description="Renders slice of scalar field with bubble edge")
parser.add_argument('files',
                    nargs='+',
                    help="List of scalar data files 'u_*.xmf'")
parser.add_argument('--force',
                    action="store_true",
                    help="Overwrite existing files")
parser.add_argument('--draft',
                    action="store_true",
                    help="Fewer lower resolution")
parser.add_argument('--res',
                    type=int,
                    nargs=2,
                    default=[1600, 800],
                    help="Image resolution in pixels")
parser.add_argument('--fieldname',
                    type=str,
                    default="",
                    help="Name of scalar field for slice color. "
                    "Defaults to prefix of filename")
parser.add_argument('--colormap',
                    type=str,
                    choices=("coolwarm", "yellow", "rainbow"),
                    default="coolwarm",
                    help="Colormap for slice color")
parser.add_argument('--levels',
                    type=int,
                    default=20,
                    help="Number of discrete values in colormap")
parser.add_argument('--bcgroup',
                    type=float,
                    default=2,
                    help="Group of boundary conditions to draw electrode")
parser.add_argument('--bubblewidth',
                    type=float,
                    default=5,
                    help="Bubble line width")
parser.add_argument('--bcwidth',
                    type=float,
                    default=8,
                    help="Electrode line width")
parser.add_argument('--bccolor',
                    type=float,
                    nargs=3,
                    default=[0.5, 0.5, 0.5],
                    help="Electrode line color")
parser.add_argument('--vmin',
                    type=float,
                    default=0,
                    help="Maximum value of scalar field")
parser.add_argument('--vmax',
                    type=float,
                    default=1,
                    help="Maximum value of scalar field")
args = parser.parse_args()

sources_ft = []
timearrays = []

if not args.fieldname:
    args.fieldname = re.findall(
        "(.*)_.*",
        os.path.splitext(os.path.basename(args.files[0]))[0])[0]

files_tu0 = args.files
source_tu0 = XDMFReader(FileNames=files_tu0)
box = paratools.GetBoundingBox(source_tu0)
boxc = (box[0] + box[1]) * 0.5
boxsize = box[1] - box[0]
source_tu0.CellArrayStatus = [args.fieldname]
(source_tu0, ), (timearray, ) = paratools.ApplyForceTime([source_tu0])
sources_ft.append(source_tu0)
timearrays.append(timearray)

files_sm = paratools.ReplaceFilename(args.files, "sm_{}.vtk")
source_sm = LegacyVTKReader(FileNames=files_sm)
(source_sm, ), (timearray, ) = paratools.ApplyForceTime([source_sm])
sources_ft.append(source_sm)
timearrays.append(timearray)

'''
files_part = paratools.ReplaceFilename(args.files, "part_{}.csv")
source_part = CSVReader(FileName=files_part)
(source_part, ), (timearray, ) = paratools.ApplyForceTime([source_part])
sources_ft.append(source_part)
timearrays.append(timearray)
'''

files_bc = paratools.ReplaceFilename(args.files[:1], "bc.vtk")
source_bc = LegacyVTKReader(FileNames=files_bc)

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.OrientationAxesVisibility = 0
def V(*v):
    return np.array(v)
renderView1.CameraPosition = boxc + V(0, 0, 10)
renderView1.CameraFocalPoint = boxc
renderView1.CameraViewUp = [0, 1, 0]
renderView1.CameraParallelScale = boxsize[1] * 0.5
renderView1.CameraParallelProjection = 1
renderView1.Background = [1] * 3
renderView1.UseLight = 0

viewsize = args.res
if args.draft:
    viewsize = [max(1, s // 2) for s in viewsize]
renderView1.ViewSize = viewsize


# Electrode from boudary conditions.
if args.bcwidth:
    slice_bc = Slice(Input=source_bc)
    slice_bc.SliceType = 'Plane'
    slice_bc.SliceType.Origin = boxc
    slice_bc.SliceType.Normal = [0.0, 0.0, 1.0]
    threshold_bc = Threshold(Input=slice_bc)
    threshold_bc.Scalars = ['CELLS', 'group']
    threshold_bc.ThresholdRange = [args.bcgroup, args.bcgroup]
    threshold_bcDisplay = Show(threshold_bc, renderView1)
    threshold_bcDisplay.Representation = 'Wireframe'
    threshold_bcDisplay.ColorArrayName = ['POINTS', '']
    threshold_bcDisplay.AmbientColor = args.bccolor
    threshold_bcDisplay.DiffuseColor = args.bccolor
    threshold_bcDisplay.LineWidth = args.bcwidth

# Concentration field.
calc_tu0 = Calculator(Input=source_tu0)
calc_tu0.AttributeType = 'Cell Data'
calc_tu0.ResultArrayName = 'u'
calc_tu0.Function = '({} - {:}) / {:}'.format(args.fieldname, args.vmin,
                                               args.vmax - args.vmin)
calc_tu0 = CellDatatoPointData(Input=calc_tu0)
calc_tu0.CellDataArraytoprocess = ['u']

slice1 = Slice(Input=calc_tu0)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = boxc
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1Display = Show(slice1, renderView1)
tu0LUT = GetColorTransferFunction('u')
tu0LUT.AutomaticRescaleRangeMode = 'Never'
tu0LUT.ScalarRangeInitialized = 1.0
if args.levels:
    tu0LUT.Discretize = 1
    tu0LUT.NumberOfTableValues = args.levels
else:
    tu0LUT.Discretize = 0

if args.colormap == "rainbow":
    tu0LUT.ColorSpace = 'RGB'
    tu0LUT.RGBPoints = [
        0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.143, 0.0, 0.0,
        0.360784313725, 0.285, 0.0, 1.0, 1.0, 0.429, 0.0, 0.501960784314, 0.0,
        0.571, 1.0, 1.0, 0.0, 0.714, 1.0, 0.380392156863, 0.0, 0.857,
        0.419607843137, 0.0, 0.0, 1.0, 0.878431372549, 0.301960784314,
        0.301960784314
    ]
elif args.colormap == "coolwarm":
    tu0LUT.ColorSpace = 'Lab'
    tu0LUT.RGBPoints = [
        0.0, 0.0, 0.0, 0.34902, 0.03125000000000003, 0.039216, 0.062745,
        0.380392, 0.06250000000000006, 0.062745, 0.117647, 0.411765,
        0.09374999999999994, 0.090196, 0.184314, 0.45098, 0.12499999999999997,
        0.12549, 0.262745, 0.501961, 0.15625, 0.160784, 0.337255, 0.541176,
        0.18750000000000003, 0.2, 0.396078, 0.568627, 0.21875000000000006,
        0.239216, 0.454902, 0.6, 0.24999999999999994, 0.286275, 0.521569,
        0.65098, 0.28124999999999994, 0.337255, 0.592157, 0.701961, 0.3125,
        0.388235, 0.654902, 0.74902, 0.34375, 0.466667, 0.737255, 0.819608,
        0.37500000000000006, 0.572549, 0.819608, 0.878431, 0.40625, 0.654902,
        0.866667, 0.909804, 0.4375, 0.752941, 0.917647, 0.941176, 0.46875,
        0.823529, 0.956863, 0.968627, 0.5, 0.988235, 0.960784, 0.901961, 0.5,
        0.941176, 0.984314, 0.988235, 0.52, 0.988235, 0.945098, 0.85098, 0.54,
        0.980392, 0.898039, 0.784314, 0.5625, 0.968627, 0.835294, 0.698039,
        0.59375, 0.94902, 0.733333, 0.588235, 0.625, 0.929412, 0.65098,
        0.509804, 0.65625, 0.909804, 0.564706, 0.435294, 0.6875, 0.878431,
        0.458824, 0.352941, 0.71875, 0.839216, 0.388235, 0.286275,
        0.7500000000000001, 0.760784, 0.294118, 0.211765, 0.78125, 0.701961,
        0.211765, 0.168627, 0.8125, 0.65098, 0.156863, 0.129412, 0.84375, 0.6,
        0.094118, 0.094118, 0.875, 0.54902, 0.066667, 0.098039,
        0.9062500000000001, 0.501961, 0.05098, 0.12549, 0.9375, 0.45098,
        0.054902, 0.172549, 0.96875, 0.4, 0.054902, 0.192157, 1.0, 0.34902,
        0.070588, 0.211765
    ]
elif args.colormap == "yellow":
    tu0LUT.ColorSpace = 'Lab'
    tu0LUT.RGBPoints = [
        0.0, 1.0, 1.0, 0.988235, 0.002, 1.0, 1.0, 0.988235,
        0.05000000000000001, 0.984314, 0.988235, 0.843137, 0.10000000000000002,
        0.988235, 0.988235, 0.741176, 0.15, 0.980392, 0.968627, 0.654902,
        0.20000000000000004, 0.980392, 0.945098, 0.576471, 0.25, 0.968627,
        0.905882, 0.486275, 0.3, 0.968627, 0.862745, 0.388235,
        0.3499999999999999, 0.960784, 0.803922, 0.286275, 0.4000000000000001,
        0.94902, 0.741176, 0.219608, 0.45, 0.941176, 0.678431, 0.14902, 0.5,
        0.929412, 0.607843, 0.094118, 0.55, 0.921569, 0.545098, 0.054902, 0.6,
        0.909804, 0.486275, 0.035294, 0.65, 0.890196, 0.411765, 0.019608,
        0.6999999999999998, 0.8, 0.305882, 0.0, 0.7500000000000001, 0.760784,
        0.239216, 0.0, 0.8000000000000002, 0.678431, 0.180392, 0.011765, 0.85,
        0.6, 0.121569, 0.023529, 0.9, 0.501961, 0.054902, 0.031373, 0.95, 0.4,
        0.039216, 0.058824, 1.0, 0.301961, 0.047059, 0.090196
    ]
else:
    assert False, "Unknown colormap=" + args.colormap
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'u']
slice1Display.LookupTable = tu0LUT

# Bubble edge.
if args.bubblewidth:
    slice_sm = Slice(Input=source_sm)
    slice_sm.SliceType = 'Plane'
    slice_sm.SliceOffsetValues = [0.0]
    slice_sm.SliceType.Origin = boxc
    slice_sm.SliceType.Normal = [0.0, 0.0, 1.0]
    slice_smDisplay = Show(slice_sm, renderView1)
    slice_smDisplay.Representation = 'Wireframe'
    slice_smDisplay.AmbientColor = [0.0, 0.0, 0.0]
    slice_smDisplay.ColorArrayName = ['POINTS', '']
    slice_smDisplay.DiffuseColor = [0.0, 0.0, 0.0]
    slice_smDisplay.LineWidth = args.bubblewidth

# FIXME: Workaround, otherwise scalar field is not shown in the first image.
paratools.SetTimeStep(1, sources_ft, timearrays)

pattern = "a_{}.png"
steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force,
                        pattern=pattern)
