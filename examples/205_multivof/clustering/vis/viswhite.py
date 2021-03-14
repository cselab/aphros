#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.1

import argparse

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help="Path to sm_*.vtk files")
parser.add_argument('--resx', type=int, default=1080, help="resolution in x")
parser.add_argument('--resy', type=int, default=1080, help="resolution in y")
parser.add_argument('--samples',
                    type=int,
                    default=100,
                    help="number of samples")
parser.add_argument('--force',
                    action="store_true",
                    help="overwrite existing files")
parser.add_argument('--camera',
                    type=str,
                    default="expclose",
                    choices=("exp", "expclose"))
parser.add_argument('--draft',
                    action="store_true",
                    help="few samples and low resolution")
parser.add_argument('--ray',
                    default=1,
                    type=int,
                    choices=(0, 1),
                    help="raytracing")
parser.add_argument('--ambient',
                    default=1,
                    type=int,
                    choices=(0, 1),
                    help="use ambient light instead of light kit")
args = parser.parse_args()

if args.draft:
    args.resx //= 2
    args.resy //= 2
    args.samples = max(1, args.samples // 10)

matpath = "tmp_materials.json"
with open(matpath, 'w') as f:
    f.write('''\
{
  "family" : "OSPRay",
  "version" : "0.0",
  "materials" : {
    "water_white": {
      "type": "Principled",
      "doubles" : {
          "baseColor" : [1, 1, 1],
          "ior" : [1.33],
          "transmissionColor" : [0.8, 0.8, 0.8],
          "transmission" : [0.95],
          "transmissionDepth" : [0.6],
          "specular" : [1],
          "diffuse" : [1]
      }
    }
  }
}
''')

if args.ray:
    materialLibrary1 = GetMaterialLibrary()
    materialLibrary1.LoadMaterials = matpath

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [args.resx, args.resy]
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0 if args.ambient else 1
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
renderView1.Background = [1] * 3

if args.camera == "exp":
    renderView1.CameraPosition = [0.711, 4.054, 1.582]
    renderView1.CameraFocalPoint = [0.711, -0.281, 0.613]
    renderView1.CameraViewUp = [0.0, 0.149, -0.988]
elif args.camera == "expclose":
    renderView1.CameraPosition = [0.711, 3.054, 1.482]
    renderView1.CameraFocalPoint = [0.711, -0.281, 0.513]
    renderView1.CameraViewUp = [0.0, 0.149, -0.988]
else:
    assert False, "Unknown camera={}".format(args.camera)

if args.ray:
    renderView1.EnableRayTracing = 1
    renderView1.BackEnd = 'OSPRay pathtracer'
    renderView1.AmbientSamples = 1
    renderView1.SamplesPerPixel = args.samples
    renderView1.OSPRayMaterialLibrary = materialLibrary1

if args.ambient and args.ray:
    light1 = CreateLight()
    light1.Coords = 'Ambient'
    light1.Intensity = 3
    renderView1.AdditionalLights = light1

if len(args.files) > 1:
    import paratools
    source_sm = LegacyVTKReader(
        FileNames=paratools.ReplaceFilename(args.files, "sm_{}.vtk"))
    sources_ft, timearrays = paratools.ApplyForceTime([source_sm])
    source_sm, = sources_ft
else:
    source_sm = LegacyVTKReader(FileNames=args.files)

surf = Calculator(Input=source_sm)
surf.ResultNormals = 1
surf.ResultArrayName = 'normals'
surf.Function = 'nn'
try:
    surf.AttributeType = 'Point Data'
except:
    pass

# Stretch boundaries to fill the view.
# Expecting bounding box between [0,0,0] and [1.5, *, 1.5]
surf = Transform(Input=surf)
surf.Transform = 'Transform'
surf.Transform.Scale = [50., 50., 50.]
surf = Calculator(Input=surf)
try:
    surf.AttributeType = 'Point Data'
except:
    pass
surf.CoordinateResults = 1
surf.ResultArrayName = 'warp'
surf.Function = '\
    (100*coordsX - 99*max(0.05, min(1.45, coordsX)))*iHat + \
    coordsY*jHat + \
    (100*coordsZ - 99*max(0.05, min(1.45, coordsZ)))*kHat'

surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = 'Surface'
surfDisplay.ColorArrayName = ['POINTS', '']
if args.ray:
    surfDisplay.OSPRayMaterial = 'water_white'
else:
    surfDisplay.Opacity = 0.75

eps = 1e-5

planecolor = [0.9] * 3
plane1 = Plane()
plane1.Origin = [-2.25, -eps, -2.25]
plane1.Point1 = [3.75, -eps, -2.25]
plane1.Point2 = [-2.25, -eps, 3.75]
plane1.XResolution = 30
plane1.YResolution = 30
plane1Display = Show(plane1, renderView1)
plane1Display.Representation = 'Surface With Edges'
plane1Display.ColorArrayName = [None, '']
plane1Display.AmbientColor = planecolor
plane1Display.DiffuseColor = planecolor
plane1Display.LineWidth = 0.4
plane1Display.EdgeColor = [0.85] * 3

if len(args.files) > 1:
    steps = paratools.GetSteps(args.files)
    paratools.SaveAnimation(steps,
                            renderView1,
                            sources_ft,
                            timearrays,
                            force=args.force)
else:
    SaveScreenshot("a.png")
