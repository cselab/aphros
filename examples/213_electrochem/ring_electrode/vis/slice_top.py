#!/usr/bin/env pvbatch

# state file generated using paraview version 5.9.0

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
import argparse
import paratools
import numpy as np

parser = argparse.ArgumentParser(
    description="Renders vorticity magnitude with LAMMPS polymers")
parser.add_argument('files', nargs='+', help="List of data files 'sm_*.vtk'")
parser.add_argument('--force',
                    action="store_true",
                    help="Overwrite existing files")
parser.add_argument('--draft',
                    action="store_true",
                    help="Fewer lower resolution")
parser.add_argument('--res',
                    type=int,
                    nargs=2,
                    default=[1920, 1080],
                    help="Image resolution in pixels")
parser.add_argument('--extent', type=float, default=1., help="Domain extent")
parser.add_argument("--plane", type=int, default=1, help="Draw bottom plane")
args = parser.parse_args()

sources_ft = []
timearrays = []

files_tu0 = paratools.ReplaceFilename(args.files, "tu0_{}.xmf")
source_tu0 = XDMFReader(FileNames=files_tu0)
source_tu0.CellArrayStatus = ['tu0']
source_tu0.GridStatus = ['mesh']
(source_tu0, ), (timearray, ) = paratools.ApplyForceTime([source_tu0])
sources_ft.append(source_tu0)
timearrays.append(timearray)

files_sm = paratools.ReplaceFilename(args.files, "sm_{}.vtk")
source_sm = LegacyVTKReader(FileNames=files_sm)
(source_sm, ), (timearray, ) = paratools.ApplyForceTime([source_sm])
sources_ft.append(source_sm)
timearrays.append(timearray)

files_part = paratools.ReplaceFilename(args.files, "part_{}.csv")
source_part = CSVReader(FileName=files_part)
(source_part, ), (timearray, ) = paratools.ApplyForceTime([source_part])
sources_ft.append(source_part)
timearrays.append(timearray)

files_bc = paratools.ReplaceFilename(args.files[:1], "bc.vtk")
source_bc = LegacyVTKReader(FileNames=files_bc)

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
ext = args.extent
def V(*v):
    return np.array(v)
renderView1.CameraPosition = V(0.5, 1.3969, 3.1080) * ext
renderView1.CameraFocalPoint = V(0.5, 0.2017, 0.5450) * ext
renderView1.CameraViewUp = V(0.0, 0.9063, -0.422) * ext
renderView1.CameraParallelScale = 0.25 * ext
renderView1.CameraParallelProjection = 1
renderView1.Background = [1] * 3

viewsize = args.res
if args.draft:
    viewsize = [max(1, s // 2) for s in viewsize]
renderView1.ViewSize = viewsize

# Bottom plane.
if args.plane:
    plane1 = Plane()
    eps = 1e-3
    plane1.Origin = [-4 * ext, -eps, -4 * ext]
    plane1.Point1 = [5 * ext, -eps, -4 * ext]
    plane1.Point2 = [-4 * ext, -eps, 5 * ext]
    plane1Display = Show(plane1, renderView1)
    plane1Display.Representation = 'Surface'
    plane1Display.ColorArrayName = [None, '']

# Boundary conditions.
threshold1 = Threshold(Input=source_bc)
threshold1.Scalars = ['CELLS', 'group']
threshold1.ThresholdRange = [2, 2]
threshold1Display = Show(threshold1, renderView1)
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['CELLS', '']
threshold1Display.DiffuseColor = [0.8, 0, 0]


# Particles.
tableToPoints2 = TableToPoints(Input=source_part)
tableToPoints2.XColumn = 'x'
tableToPoints2.YColumn = 'y'
tableToPoints2.ZColumn = 'z'
glyph2 = Glyph(Input=tableToPoints2, GlyphType='Sphere')
glyph2.OrientationArray = ['POINTS', 'No orientation array']
glyph2.ScaleArray = ['POINTS', 'r']
glyph2.GlyphTransform = 'Transform2'
glyph2.GlyphMode = 'All Points'
glyph2.GlyphType.Radius = 1.0
glyph2.GlyphType.ThetaResolution = 20
glyph2.GlyphType.PhiResolution = 20
glyph2Display = Show(glyph2, renderView1)
glyph2Display.Representation = 'Surface'
glyph2Display.AmbientColor = [0.0, 0.0, 0.0]
glyph2Display.ColorArrayName = ['POINTS', '']
glyph2.ScaleFactor = 1

# Concentration field.
slice1 = Slice(Input=source_tu0)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = V(0.5, 0.5, 0.5) * ext
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1 = CellDatatoPointData(Input=slice1)
slice1.CellDataArraytoprocess = ['tu0']
slice1Display = Show(slice1, renderView1)
tu0LUT = GetColorTransferFunction('tu0')
tu0LUT.AutomaticRescaleRangeMode = 'Never'
tu0LUT.RGBPoints = [
    0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.08579999999999999,
    0.0, 0.0, 0.360784313725, 0.17099999999999999, 0.0, 1.0, 1.0,
    0.25739999999999996, 0.0, 0.501960784314, 0.0, 0.34259999999999996, 1.0,
    1.0, 0.0, 0.42839999999999995, 1.0, 0.380392156863, 0.0, 0.5142,
    0.419607843137, 0.0, 0.0, 0.6, 0.878431372549, 0.301960784314,
    0.301960784314
]
tu0LUT.ColorSpace = 'RGB'
tu0LUT.ScalarRangeInitialized = 1.0
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'tu0']
slice1Display.LookupTable = tu0LUT
slice1Display.Ambient = 0.3

contour1 = Contour(Input=slice1)
contour1.ContourBy = ['POINTS', 'tu0']
contour1.Isosurfaces = np.linspace(0, 0.6, 21)
contour1Display = Show(contour1, renderView1)
contour1Display.Representation = 'Wireframe'
contour1Display.LineWidth = 2
contour1Display.ColorArrayName = ['POINTS', '']
contour1Display.DiffuseColor = [1, 1, 1]

# Bubble surfaces.
calcNormals = Calculator(Input=source_sm)
calcNormals.ResultNormals = 1
calcNormals.AttributeType = 'Point Data'
calcNormals.ResultArrayName = 'normals'
calcNormals.Function = 'nn'
calcNormalsDisplay = Show(calcNormals, renderView1)
calcNormalsDisplay.Representation = 'Surface'
if hasattr(calcNormalsDisplay, "SelectNormalArray"):
    calcNormalsDisplay.SelectNormalArray = 'normals'
calcNormalsDisplay.AmbientColor = [1, 1, 1]
calcNormalsDisplay.ColorArrayName = ['POINTS', '']
calcNormalsDisplay.DiffuseColor = [1, 1, 1]
calcNormalsDisplay.ColorArrayName = [None, '']

# FIXME: Workaround, otherwise `tu0` not shown in the first image.
paratools.SetTimeStep(1, sources_ft, timearrays)

pattern = "a_{}.png"
steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force,
                        pattern=pattern)

