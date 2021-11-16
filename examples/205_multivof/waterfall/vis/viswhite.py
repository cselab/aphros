#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.1

import argparse

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

parser = argparse.ArgumentParser()
parser.add_argument(
    'files',
    nargs='+',
    default='/scratch/snx3000/karnakov/pub/waterfall/sm_0282.vtk',
    help="Path to a sm_*.vtk file")
parser.add_argument('--resx', type=int, default=1920, help="resolution in x")
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
                    default="aps",
                    choices=("aps", "close", "close2", "cover"))
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
          "transmissionColor" : [0.42, 0.74, 0.95],
          "transmission" : [0.9],
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

if args.camera == "aps":
    renderView1.CameraPosition = [
        1.645754077535257, 1.1106446111614987, 2.977383845728232
    ]
    renderView1.CameraFocalPoint = [
        0.48999991817830235, -0.26932721490693573, -1.3056414919409265
    ]
    renderView1.CameraViewUp = [
        -0.0893981734829278, 0.9547941520217919, -0.2835067791833273
    ]
elif args.camera == "close":
    renderView1.CameraPosition = [
        1.6866822392533614, 1.571758849156841, 2.0958321780122997
    ]
    renderView1.CameraFocalPoint = [
        0.36452963249954384, -0.9285898276110117, -1.5898923016211817
    ]
    renderView1.CameraViewUp = [
        -0.2803861781098946, 0.838087397250151, -0.46796699210029125
    ]
elif args.camera == "close2":
    renderView1.CameraPosition = [
        2.0838886233689013, 0.9163047538734606, 1.6684874774246188
    ]
    renderView1.CameraFocalPoint = [
        -1.0060111943225836, -0.2918065772628064, -1.583801268818597
    ]
    renderView1.CameraViewUp = [
        -0.1847806328745186, 0.9655687359221463, -0.18312054479004014
    ]
elif args.camera == "cover":
    renderView1.CameraPosition = [1.331, 3.0769, 3.8146]
    renderView1.CameraFocalPoint = [0.95301, 0.71242, -0.1667]
    renderView1.CameraViewUp = [-0.11556, 0.85882, -0.49908]
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
    light1.Intensity = 4.5
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

surfDisplay = Show(surf, renderView1)
surfDisplay.Representation = 'Surface'
surfDisplay.ColorArrayName = ['POINTS', '']
if args.ray:
    surfDisplay.OSPRayMaterial = 'water_white'

eps = 1e-5

# bottom plane
plane1color = [0.5] * 3
width = 100
plane1 = Plane()
plane1.Origin = [-width, -eps, -width]
plane1.Point1 = [width, -eps, -width]
plane1.Point2 = [-width, -eps, width]
plane1Display = Show(plane1, renderView1)
plane1Display.Representation = 'Surface'
plane1Display.AmbientColor = plane1color
plane1Display.DiffuseColor = plane1color
plane1Display.ColorArrayName = [None, '']

# back plane
plane2color = [0.5] * 3
plane2 = Plane()
plane2.Origin = [-width, 0, -eps]
plane2.Point1 = [width, 0, -eps]
plane2.Point2 = [-width, width, -eps]
plane2Display = Show(plane2, renderView1)
plane2Display.Representation = 'Surface'
plane2Display.AmbientColor = plane2color
plane2Display.DiffuseColor = plane2color
plane2Display.ColorArrayName = [None, '']

if len(args.files) > 1:
    steps = paratools.GetSteps(args.files)
    paratools.SaveAnimation(steps,
                            renderView1,
                            sources_ft,
                            timearrays,
                            force=args.force)
else:
    SaveScreenshot("a.png")
