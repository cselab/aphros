#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np

import paratools

parser = argparse.ArgumentParser(
    description=
    "Renders slice with gas concentration and isolines of electric field.")
parser.add_argument('files',
                    nargs='+',
                    help="List of data files 'elpot_*.xmf'")
parser.add_argument('--force',
                    action="store_true",
                    help="Overwrite existing files")
parser.add_argument('--draft', action="store_true", help="Lower resolution")
parser.add_argument('--ebvf',
                    type=int,
                    default=1,
                    help="Read embed volume fraction from ebvf_0000.xmf"
                    " and apply threshold to hide excluded cells")
parser.add_argument('--resy',
                    type=int,
                    default="1000",
                    help="Output image height")
args = parser.parse_args()

if args.draft:
    args.resy //= 2

source_bcvtk = LegacyVTKReader(
    FileNames=paratools.ReplaceFilename(args.files[:1], "bc.vtk"))

source_elpot = XDMFReader(
    FileNames=paratools.ReplaceFilename(args.files, "elpot_{}.xmf"))
source_elpot.CellArrayStatus = ['elpot']
source_elpot.GridStatus = ['Grid_10220']

source_tu0 = XDMFReader(
    FileNames=paratools.ReplaceFilename(args.files, "tu0_{}.xmf"))
source_tu0.CellArrayStatus = ['tu0']
source_tu0.GridStatus = ['Grid_10220']

source_tu1 = XDMFReader(
    FileNames=paratools.ReplaceFilename(args.files, "tu1_{}.xmf"))
source_tu1.CellArrayStatus = ['tu1']
source_tu1.GridStatus = ['Grid_10220']

sources_ft, timearrays = paratools.ApplyForceTime(
    [source_elpot, source_tu0, source_tu1])
source_elpot, source_tu0, source_tu1 = sources_ft

if args.ebvf:
    source_ebvf = XDMFReader(
        FileNames=paratools.ReplaceFilename(args.files, "ebvf_{}.xmf"))
    source_ebvf.CellArrayStatus = ['ebvf']
    source_ebvf.GridStatus = ['Grid_10220']
    (source_ebvf, ), (timearray, ) = paratools.ApplyForceTime([source_ebvf])
    sources_ft.append(source_ebvf)
    timearrays.append(timearray)

box = paratools.GetBoundingBox(source_bcvtk)
boxc = (box[0] + box[1]) * 0.5
boxsize = box[1] - box[0]


def ceil(a, d):
    return (int(a) + d - 1) // d * d


div = 4
viewsize = boxsize * args.resy / boxsize[1]
viewsize = [ceil(viewsize[0], div), ceil(viewsize[1], div)]

renderView1 = CreateView('RenderView')
renderView1.ViewSize = viewsize
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.CameraPosition = [boxc[0], boxc[1], 10]
renderView1.CameraFocalPoint = [boxc[0], boxc[1], 0]
renderView1.CameraParallelScale = boxsize[1] * 0.5 * 1.05
renderView1.CameraParallelProjection = 1
renderView1.Background = [1] * 3

slice_bc = Slice(Input=source_bcvtk)
slice_bc.SliceType = 'Plane'
slice_bc.SliceOffsetValues = [0.0]
slice_bc.SliceType.Origin = boxc
slice_bc.SliceType.Normal = [0.0, 0.0, 1.0]
slice_bcDisplay = Show(slice_bc, renderView1)
slice_bcDisplay.Representation = 'Wireframe'
slice_bcDisplay.ColorArrayName = ['POINTS', '']
slice_bcDisplay.LineWidth = 5.0
slice_bcDisplay.DiffuseColor = [0.0, 0.0, 0.0]
slice_bcDisplay.AmbientColor = [0.0, 0.0, 0.0]

attr = [source_tu0, source_tu1, source_elpot]
if args.ebvf:
    attr.append(source_ebvf)
append_attr = AppendAttributes(Input=attr)
append_attr = CellDatatoPointData(Input=append_attr)
attr_names = ['elpot', 'tu0', 'tu1']
if args.ebvf:
    attr_names.append('ebvf')
append_attr.CellDataArraytoprocess = attr_names
slice_append = Slice(Input=append_attr)
slice_append.SliceType = 'Plane'
slice_append.SliceOffsetValues = [0.0]
slice_append.SliceType.Origin = boxc
slice_append.SliceType.Normal = [0.0, 0.0, 1.0]
if args.ebvf:
    slice_append = IsoVolume(Input=slice_append)
    slice_append.InputScalars = ['POINTS', 'ebvf']
    slice_append.ThresholdRange = [0.5, 2.0]

contour_elpot = Contour(Input=slice_append)
contour_elpot.ContourBy = ['POINTS', 'elpot']
contour_elpot.Isosurfaces = np.linspace(0.05, 0.95, 31)
contour_elpot.PointMergeMethod = 'Uniform Binning'
contour_elpotDisplay = Show(contour_elpot, renderView1)
contour_elpotDisplay.Representation = 'Surface'
contour_elpotDisplay.ColorArrayName = ['POINTS', '']
contour_elpotDisplay.DiffuseColor = [0.0, 0.0, 0.0]
contour_elpotDisplay.LineWidth = 3.0

calc_tu = Calculator(Input=slice_append)
calc_tu.AttributeType = 'Point Data'
calc_tu.ResultArrayName = 'tu'
calc_tu.Function = 'tu0+tu1'
calc_tuDisplay = Show(calc_tu, renderView1)
tuLUT = GetColorTransferFunction('tu')
tuLUT.AutomaticRescaleRangeMode = 'Never'
tuLUT.RGBPoints = [
    0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.0715, 0.0, 0.0,
    0.360784313725, 0.1425, 0.0, 1.0, 1.0, 0.2145, 0.0, 0.501960784314, 0.0,
    0.2855, 1.0, 1.0, 0.0, 0.357, 1.0, 0.380392156863, 0.0, 0.4285,
    0.419607843137, 0.0, 0.0, 0.5, 0.878431372549, 0.301960784314,
    0.301960784314
]
tuLUT.ColorSpace = 'RGB'
tuLUT.ScalarRangeInitialized = 1.0
calc_tuDisplay.Representation = 'Surface'
calc_tuDisplay.ColorArrayName = ['POINTS', 'tu']
calc_tuDisplay.LookupTable = tuLUT

threshold_electrodes = Threshold(registrationName='threshold_electrodes',
                                 Input=slice_bc)
threshold_electrodes.Scalars = ['CELLS', 'group']
threshold_electrodes.ThresholdRange = [3.0, 4.0]
threshold_electrodesDisplay = Show(threshold_electrodes, renderView1)
threshold_electrodesDisplay.Representation = 'Wireframe'
threshold_electrodesDisplay.AmbientColor = [1.0, 0.0, 1.0]
threshold_electrodesDisplay.ColorArrayName = ['POINTS', '']
threshold_electrodesDisplay.DiffuseColor = [1.0, 0.0, 1.0]
threshold_electrodesDisplay.LineWidth = 8.0
threshold_electrodesDisplay.Position = [0.0, 0.0, 0.1]
threshold_electrodesDisplay.PolarAxes.Translation = [0.0, 0.0, 0.1]

steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)
