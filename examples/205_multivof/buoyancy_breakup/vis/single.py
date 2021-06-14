#!/usr/bin/env pvbatch

import argparse
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import os


parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help="Path to sm_*.vtk files")
parser.add_argument('outdir', default='.', help="Path to output directory")
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

surf = LegacyVTKReader(FileNames=args.files)

#d3 category10
vhex = [
    "1f77b4",
    "ff7f0e",
    "2ca02c",
    "d62728",
    "9467bd",
    "8c564b",
    "e377c2",
    "7f7f7f",
    "bcbd22",
    "17becf",
]


def rgb(h):
    return list(int(h[i:i + 2], 16) / 255. for i in (0, 2, 4))

lw = 8

surfDisplay = Show(surf, renderView1, 'GeometryRepresentation')
surfDisplay.Representation = 'Wireframe'
surfDisplay.AmbientColor = rgb(vhex[0])
surfDisplay.ColorArrayName = ['POINTS', '']
surfDisplay.DiffuseColor = rgb(vhex[0])
surfDisplay.LineWidth = lw
surfDisplay.Position = [0.0, 0.0, 0.0]

tk = GetTimeKeeper()
for i, f in enumerate(args.files):
    tk.Time = i
    path = os.path.join(args.outdir,
                        os.path.splitext(os.path.basename(f))[0] + '.png')
    print(path)
    SaveScreenshot(path)
