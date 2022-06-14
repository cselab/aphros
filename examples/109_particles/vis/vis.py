#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np

import paratools

parser = argparse.ArgumentParser(description="Renders particles.")
parser.add_argument('files',
                    nargs='+',
                    help="list of data files 'particles_*.csv'")
parser.add_argument('--force',
                    action="store_true",
                    help="overwrite existing files")
parser.add_argument('--omz_max',
                    type=float,
                    default=10,
                    help="maximum vorticity omz, range will be [-max,max]")
parser.add_argument('--resy',
                    type=int,
                    nargs=1,
                    default=1080,
                    help="output image resolution in y")
parser.add_argument('--colormap',
                    type=str,
                    default='geo_d2',
                    choices=['geo_d2', 'geo_br'],
                    help="maximum vorticity omz, range will be [-max,max]")
args = parser.parse_args()

sources_ft = []
timearrays = []

source_csv = CSVReader(FileName=args.files)
(source_csv, ), (timearray, ) = paratools.ApplyForceTime([source_csv])
sources_ft.append(source_csv)
timearrays.append(timearray)

renderView1 = CreateView('RenderView')
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.ViewSize = [args.resy, args.resy]
renderView1.CameraPosition = [0.5, 0.5, 10]
renderView1.CameraFocalPoint = [0.5, 0.5, 0]
renderView1.CameraParallelScale = 0.5
renderView1.CameraParallelProjection = 1
renderView1.Background = [1] * 3

tableToPoints1 = TableToPoints(Input=source_csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'
tableToPoints1Display = Show(tableToPoints1, renderView1)
tableToPoints1Display.Representation = 'Points'
tableToPoints1Display.AmbientColor = [0.0, 0.0, 0.0]
tableToPoints1Display.ColorArrayName = [None, '']
tableToPoints1Display.DiffuseColor = [0.0, 0.0, 0.0]
tableToPoints1Display.PointSize = 3.

steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)
