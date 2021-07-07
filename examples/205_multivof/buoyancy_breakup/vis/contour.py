#!/usr/bin/env pvbatch

import argparse
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import os
import sys

def printerr(m):
    sys.stderr.write('{:}\n'.format(m))
    sys.stderr.flush()


parser = argparse.ArgumentParser()
parser.add_argument('--files0', nargs='*', help="Paths to sm_*.vtk files to render with color C0")
parser.add_argument('--files1', nargs='*', help="Paths to sm_*.vtk files to render with color C1")
parser.add_argument('--files2', nargs='*', help="Paths to sm_*.vtk files to render with color C2")
parser.add_argument('--files3', nargs='*', help="Paths to sm_*.vtk files to render with color C3")
parser.add_argument('--lw', default=4, help="Line width")
parser.add_argument('--force', action='store_true', help="Force overwrite")
parser.add_argument('--outdir', default='.', help="Path to output directory")
args = parser.parse_args()

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1080, 1080]
renderView1.InteractionMode = '2D'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 0.5, 0]
renderView1.CameraPosition = [0.5, 0.5, 3.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.5
renderView1.CameraParallelProjection = 1
renderView1.Background = [1., 1., 1.]
renderView1.UseLight = 0

# https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit/tree/master/Colours
colorscheme = [
    "#FF1F5B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E", "#F28522",
    "#A0B1BA", "#A6761D", "#E9002D", "#FFAA00", "#00B000"
]

def rgb(h):
    return list(int(h[1:][i:i + 2], 16) / 255. for i in (0, 2, 4))

for i in range(4):
    files = eval("args.files" + str(i))
    if not files:
        continue
    surf = LegacyVTKReader(FileNames=files)
    surfDisplay = Show(surf, renderView1)
    surfDisplay.Representation = 'Wireframe'
    surfDisplay.AmbientColor = rgb(colorscheme[i])
    surfDisplay.ColorArrayName = ['POINTS', '']
    surfDisplay.DiffuseColor = rgb(colorscheme[i])
    surfDisplay.LineWidth = args.lw

tk = GetTimeKeeper()
for i, f in enumerate(args.files0):
    path = os.path.join(args.outdir,
                        os.path.splitext(os.path.basename(f))[0] + '.png')
    if not args.force and os.path.isfile(path):
        printerr("skip existing '{}'".format(path))
        continue
    tk.Time = i
    printerr(path)
    SaveScreenshot(path)

