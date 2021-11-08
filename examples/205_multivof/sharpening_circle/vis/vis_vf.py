#!/usr/bin/env pvbatch

# state file generated using paraview version 5.9.0

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import argparse
import os
import re
import paratools


parser = argparse.ArgumentParser(
    description="Renders interface shapes from s_*.vtk files.")
parser.add_argument('files', nargs='*', help="list of data files 's_*.vtk'")
parser.add_argument('--force',
                    action="store_true",
                    help="overwrite existing files")
args = parser.parse_args()

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1080, 1080]
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.CameraPosition = [2, 2, 4]
renderView1.CameraFocalPoint = [2, 2, 0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelProjection = 1
renderView1.CameraParallelScale = 1.5
renderView1.Background = [1] * 3

# https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit/tree/master/Colours
clhex_geo = [
    "FF1F5B", "00CD6C", "009ADE", "AF58BA", "FFC61E", "F28522", "A0B1BA",
    "A6761D", "E9002D", "FFAA00", "00B000"
]


def rgb(h):
    return list(int(h[i:i + 2], 16) / 255. for i in (0, 2, 4))


steps = paratools.ReplaceFilename(args.files, '{}', keep_dir=False)
source_s = LegacyVTKReader(
    FileNames=paratools.ReplaceFilename(args.files, 's_{}.vtk'))
source_vf = XDMFReader(
    FileNames=paratools.ReplaceFilename(args.files, 'vf_{}.xmf'))
source_vf.CellArrayStatus = ['vf']
sources_ft, timearrays = paratools.ApplyForceTime([source_s, source_vf])
source_s, source_vf = sources_ft

vfDisplay = Show(source_vf, renderView1)
vfLUT = GetColorTransferFunction('vf')
vfLUT.AutomaticRescaleRangeMode = 'Never'
vfLUT.RGBPoints = [0.0] + [1] * 3 + [1.0] + rgb(clhex_geo[2])
vfLUT.ColorSpace = 'RGB'
vfLUT.NanColor = [1.0, 0.0, 0.0]
vfLUT.Discretize = 0
vfLUT.ScalarRangeInitialized = 1.0
vfPWF = GetOpacityTransferFunction('vf')
vfPWF.ScalarRangeInitialized = 1
vfDisplay.Representation = 'Surface'
vfDisplay.ColorArrayName = ['CELLS', 'vf']
vfDisplay.LookupTable = vfLUT

surfDisplay = Show(source_s, renderView1)
surfDisplay.Representation = 'Wireframe'
surfDisplay.AmbientColor = [0.0, 0.0, 0.0]
surfDisplay.ColorArrayName = ['POINTS', '']
surfDisplay.DiffuseColor = [0.0, 0.0, 0.0]
surfDisplay.LineWidth = 8

paratools.SetTimeStep(1, sources_ft, timearrays)

paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)
