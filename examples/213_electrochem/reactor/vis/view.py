#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np

import paratools

parser = argparse.ArgumentParser(description="3D view of reactor")
parser.add_argument('files', nargs='+', help="List of data files 'tu0_*.xmf'")
parser.add_argument('--force',
                    action="store_true",
                    help="Overwrite existing files")
parser.add_argument('--fieldname',
                    type=str,
                    default="",
                    help="Name of scalar field for slice color. "
                    "Defaults to prefix of filename")
parser.add_argument('--part',
                    type=int,
                    default=1,
                    help="Particles from part_*.csv")
parser.add_argument('--volume', type=int, default=0, help="Volume rendering")
parser.add_argument('--colormap',
                    type=str,
                    choices=("coolwarm", "geo", "yellow", "rainbow"),
                    default="geo",
                    help="Colormap for slice color")
parser.add_argument('--camera',
                    type=str,
                    choices=("closeup", "default"),
                    default="default",
                    help="Camera position")
parser.add_argument('--vmin',
                    type=float,
                    default=0,
                    help="Maximum value of scalar field")
parser.add_argument('--vmax',
                    type=float,
                    default=1,
                    help="Maximum value of scalar field")
parser.add_argument('--draft', action="store_true", help="Lower resolution")
parser.add_argument('--resy',
                    type=int,
                    default="1080",
                    help="Output image height")
args = parser.parse_args()

if args.draft:
    args.resy //= 2

source_stepwise = XDMFReader(
    FileNames=paratools.ReplaceFilename(args.files[:1], "stepwise.xmf"))
source_stepwise.CellArrayStatus = ['stepwise']

source_vm = XDMFReader(
    FileNames=paratools.ReplaceFilename(args.files, "vm_{}.xmf"))

source_tu0 = XDMFReader(
    FileNames=paratools.ReplaceFilename(args.files, "tu0_{}.xmf"))

source_tu1 = XDMFReader(
    FileNames=paratools.ReplaceFilename(args.files, "tu1_{}.xmf"))

source_surf = LegacyVTKReader(
    FileNames=paratools.ReplaceFilename(args.files, "sm_{}.vtk"))

sources_ft, timearrays = paratools.ApplyForceTime(
    [source_vm, source_tu0, source_tu1, source_surf])
source_vm, source_tu0, source_tu1, source_surf = sources_ft

if args.part:
    source_part = CSVReader(
        FileName=paratools.ReplaceFilename(args.files, "part_{}.csv"))
    (source_part, ), (timearray, ) = paratools.ApplyForceTime([source_part])
    sources_ft.append(source_part)
    timearrays.append(timearray)

attr = [source_tu0, source_tu1, source_vm, source_stepwise]
append_attr = AppendAttributes(Input=attr)
append_attr = Threshold(Input=append_attr)
append_attr.Scalars = ['CELLS', 'stepwise']
append_attr.ThresholdRange = [0.5, 1.0]

paratools.SetTimeStep(1, sources_ft, timearrays)
paratools.SetTimeStep(0, sources_ft, timearrays)
box = paratools.GetBoundingBox(append_attr)
boxc = (box[0] + box[1]) * 0.5
boxsize = box[1] - box[0]
print("Bounding box: {:} {:}".format(list(box[0]), list(box[1])))


def V(*v):
    return np.array(v)


renderView1 = CreateView('RenderView')
renderView1.OrientationAxesVisibility = 0
renderView1.Background = [1] * 3
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5

boxcref = V(4, 1, 0.226) / 2

if args.camera == "closeup":
    renderView1.CameraPosition = (V(-0.069693, 4.507, 6.3043) -
                                  boxcref) * boxsize[1] + boxc
    renderView1.CameraFocalPoint = (V(2.2923, 0.5199, -0.18512) -
                                    boxcref) * boxsize[1] + boxc
    renderView1.CameraViewUp = V(0.17101, 0.86603, -0.46985) * boxsize[1]
    renderView1.CameraParallelScale = 0.31854 * boxsize[1]
    renderView1.CameraParallelProjection = 1
    renderView1.ViewSize = [1920 * args.resy // 1080, args.resy]
else:
    renderView1.CameraPosition = (V(-1.9269, 1.8848, 6.9265) -
                                  boxcref) * boxsize[1] + boxc
    renderView1.CameraFocalPoint = (V(2, 0.5, 0.125) -
                                    boxcref) * boxsize[1] + boxc
    renderView1.CameraViewUp = V(0.086824, 0.98481, -0.15038) * boxsize[1]
    renderView1.CameraParallelScale = 1.0258 * boxsize[1]
    renderView1.CameraParallelProjection = 1
    renderView1.ViewSize = [1920 * args.resy // 1080, args.resy]

# Concentration field.
calc_tu = Calculator(Input=append_attr)
calc_tu.AttributeType = 'Cell Data'
calc_tu.ResultArrayName = 'tu'
calc_tu.Function = 'tu0+tu1'
append_attr = calc_tu

calc_u = Calculator(Input=append_attr)
calc_u.AttributeType = 'Cell Data'
calc_u.ResultArrayName = 'u'
calc_u.Function = '({} - {:}) / {:}'.format(args.fieldname, args.vmin,
                                            args.vmax - args.vmin)
append_attr = calc_u

tu0LUT = GetColorTransferFunction('u')
tu0LUT.AutomaticRescaleRangeMode = 'Never'
tu0LUT.ScalarRangeInitialized = 1.0
if args.colormap == "rainbow":
    tu0LUT.ColorSpace = 'RGB'
    tu0LUT.RGBPoints = [
        0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.143, 0.0, 0.0,
        0.360784313725, 0.285, 0.0, 1.0, 1.0, 0.429, 0.0, 0.501960784314, 0.0,
        0.571, 1.0, 1.0, 0.0, 0.714, 1.0, 0.380392156863, 0.0, 0.857,
        0.419607843137, 0.0, 0.0, 1.0, 0.878431372549, 0.301960784314,
        0.301960784314
    ]
elif args.colormap == "geo":
    tu0LUT.ColorSpace = 'RGB'
    tu0LUT.RGBPoints = [
        0.0,
        0.0,
        0.6039215686274509,
        0.8705882352941177,  #
        0.5,
        1.0,
        1.0,
        1.0,  #
        1.0,
        1.0,
        0.12156862745098039,
        0.3568627450980392
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

if args.volume:
    tu0PWF = GetOpacityTransferFunction('u')
    tu0PWF.Points = [0, 0, 0.5, 0, 1, 1, 0.5, 0]
    tu0PWF.ScalarRangeInitialized = 1
    volume_tu = CellDatatoPointData(Input=append_attr)
    volume_tu.CellDataArraytoprocess = ['u']
    tu0Display = Show(volume_tu, renderView1)
    tu0Display.Representation = 'Volume'
    tu0Display.ColorArrayName = ['POINTS', 'u']
    tu0Display.LookupTable = tu0LUT
    tu0Display.ScalarOpacityUnitDistance = 0.005
    tu0Display.ScalarOpacityFunction = tu0PWF
    tu0Display.OpacityArrayName = ['POINTS', 'u']
else:  # Slice.
    slice_tu = Slice(Input=append_attr)
    slice_tu.SliceType = 'Plane'
    slice_tu.SliceType.Origin = boxc
    slice_tu.SliceType.Normal = [0, 0, 1]
    slice_tu = CellDatatoPointData(Input=slice_tu)
    slice_tuDisplay = Show(slice_tu, renderView1)
    slice_tuDisplay.Opacity = 0.6
    slice_tuDisplay.Ambient = 0.3
    slice_tuDisplay.Representation = 'Surface'
    slice_tuDisplay.ColorArrayName = ['POINTS', 'u']
    slice_tuDisplay.LookupTable = tu0LUT

extract_bc = ExtractSurface(Input=append_attr)
clip_bc = Clip(Input=extract_bc)
clip_bc.ClipType = 'Plane'
clip_bc.ClipType.Origin = boxc + [0, 0, boxsize[2] * 0.49]
clip_bc.ClipType.Normal = [0, 0, 1]
clip_bcDisplay = Show(clip_bc, renderView1)
clip_bcDisplay.Representation = 'Surface'
clip_bcDisplay.ColorArrayName = ['POINTS', '']

calcNormals = Calculator(Input=source_surf)
calcNormals.ResultNormals = 1
calcNormals.AttributeType = 'Point Data'
calcNormals.ResultArrayName = 'normals'
calcNormals.Function = 'nn'
calcNormalsDisplay = Show(calcNormals, renderView1)
calcNormalsDisplay.Representation = 'Surface'
if hasattr(calcNormalsDisplay, "SelectNormalArray"):
    calcNormalsDisplay.SelectNormalArray = 'normals'
calcNormalsDisplay.AmbientColor = [1, 0.66, 0]
calcNormalsDisplay.ColorArrayName = ['POINTS', '']
calcNormalsDisplay.DiffuseColor = [1, 0.66, 0]
calcNormalsDisplay.ColorArrayName = [None, '']

if args.part:
    part_points = TableToPoints(Input=source_part)
    part_points.XColumn = 'x'
    part_points.YColumn = 'y'
    part_points.ZColumn = 'z'

    part_glyph = Glyph(Input=part_points, GlyphType='Sphere')
    part_glyph.OrientationArray = ['POINTS', 'No orientation array']
    part_glyph.GlyphTransform = 'Transform2'
    part_glyph.GlyphMode = 'All Points'
    part_glyph.ScaleArray = ['POINTS', 'r']
    part_glyph.ScaleFactor = 1
    part_glyph.GlyphType.Radius = 1
    part_glyph.GlyphType.ThetaResolution = 24
    part_glyph.GlyphType.PhiResolution = 24
    partDisplay = Show(part_glyph, renderView1)
    partDisplay.Representation = 'Surface'
    partDisplay.DiffuseColor = [0, 0.66, 0]
    partDisplay.AmbientColor = [0, 0.66, 0]
    partDisplay.ColorArrayName = [None, '']
    partDisplay.SelectOrientationVectors = 'None'
    partDisplay.SelectScaleArray = 'None'
    partDisplay.SetScaleArray = [None, '']
    partDisplay.OpacityArray = [None, '']

steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)
