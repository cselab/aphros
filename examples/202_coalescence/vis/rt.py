#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0


from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import argparse
import os
import re

def Log(s):
    sys.stderr.write(str(s) + '\n')

parser = argparse.ArgumentParser(
    description="Renders interface shapes from sm_*.vtk files.")
parser.add_argument('files',
                    nargs='*',
                    help="list of data files 'sm_*.vtk'")
parser.add_argument('--force',
                    action="store_true",
                    help="overwrite existing files")
parser.add_argument('--draft',
                    action="store_true",
                    help="render at reduced resolution (half) and fewer samples")
parser.add_argument('--res',
                    type=int,
                    nargs=2,
                    default=[1920, 1080],
                    help="image resolution")
parser.add_argument('--camera',
                    type=str,
                    default="std",
                    choices=["std", "vertical"])
parser.add_argument('--preset',
                    type=str,
                    default="std",
                    choices=["", "std", "vertical"])
parser.add_argument('--samples',
                    default=50,
                    type=int,
                    help="number of samples for pathtracer")
args = parser.parse_args()

if args.preset == "std":
    args.camera = "std"
    args.res = [1920, 1080]
elif args.preset == "vertical":
    args.camera = "vertical"
    args.res = [960, 1080]

light1 = CreateLight()
light1.Coords = 'Ambient'
light1.Intensity = 1

materials_json = "materials.json"
with open(materials_json, 'w') as f:
    f.write('''{
  "family" : "OSPRay",
  "version" : "0.0",
  "materials" : {
    "bubble" : {
      "type" : "Alloy",
      "doubles" : {
        "color" : BUBBLE_COLOR,
        "edgeColor" : BUBBLE_COLOR,
        "roughness" : [0.2]
      }
    },
    "mirror" : {
      "type" : "Alloy",
      "doubles" : {
        "color" : MIRROR_COLOR,
        "edgeColor" : MIRROR_COLOR,
        "roughness" : [0.2]
      }
    },
    "scratched": {
      "type": "Principled",
      "doubles" : {
          "roughnessMap.scale" : [0.3, 0.3],
          "baseColorMap.scale" : [0.3, 0.3],
          "normalMap.scale" : [0.3, 0.3],
          "metallic" : [1.0]
      },
      "textures" : {
          "normalMap" : "TEXTURE1",
          "baseColorMap" : "TEXTURE2",
          "roughnessMap" : "TEXTURE3"
      }
    }
  }
}
'''#
            .replace("BUBBLE_COLOR", str([0.3] * 3)) #
            .replace("MIRROR_COLOR", str([0.7] * 3)) #
            .replace("TEXTURE1", os.path.abspath("metal_Normal.jpg"))
            .replace("TEXTURE2", os.path.abspath("metal_Base_Color.jpg"))
            .replace("TEXTURE3", os.path.abspath("metal_Roughness.jpg"))
    )
# XXX textures are part of ParaView distrubution
#    copy or symlink them from 'materials/'
materialLibrary = GetMaterialLibrary()
materialLibrary.LoadMaterials = materials_json

if args.draft:
    k = 2
    args.res[0] //= k
    args.res[1] //= k
    args.samples //= 10

renderView = CreateView('RenderView')
renderView.ViewSize = args.res
renderView.OrientationAxesVisibility = 0
renderView.UseLight = 1
renderView.KeyLightWarmth = 0.5
#renderView.KeyLightIntensity = 1
renderView.FillLightWarmth = 0.5

if args.camera == "std":
    renderView.CameraPosition = [4, 2.55, 14.4]
    renderView.CameraFocalPoint = [4, 1.89, 5]
    renderView.CameraViewUp = [0, 1, -0.07]
    renderView.CameraParallelProjection = 1
    renderView.CameraParallelScale = 2.45
elif args.camera == "vertical":
    renderView.CameraPosition = [4, 4.6778164394271835, 13.88409281728644]
    renderView.CameraFocalPoint = [4, 2.676401642283468, 4.672994141787723]
    renderView.CameraViewUp = [0, 0.977198335864855, -0.21232854820527072]
    renderView.CameraParallelProjection  = 1
    renderView.CameraParallelScale = 2.701125
else:
    assert False, "Unknown camera=" + args.camera


renderView.Background = [0.3] * 3
renderView.EnableRayTracing = 1
renderView.BackEnd = 'OSPRay pathtracer'
renderView.AmbientSamples = 0
renderView.SamplesPerPixel = args.samples
renderView.AdditionalLights = [light1]
renderView.OSPRayMaterialLibrary = materialLibrary

sm = LegacyVTKReader(FileNames=args.files)
surf = Calculator(Input=sm)
surf.ResultNormals = 1
surf.AttributeType = 'Point Data'
surf.ResultArrayName = 'normals'
surf.Function = 'nn'
surfDisplay = Show(surf, renderView, 'GeometryRepresentation')
surfDisplay.Representation = 'Surface'
surfDisplay.ColorArrayName = ['POINTS', '']
surfDisplay.OSPRayScaleArray = 'Normals'
surfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
surfDisplay.OSPRayMaterial = 'scratched'


bottom = Plane()
bottom.Origin = [-10.0, 0.0, -10.0]
bottom.Point1 = [18.0, 0.0, -10.0]
bottom.Point2 = [-10.0, 0.0, 14.0]
bottomDisplay = Show(bottom, renderView, 'GeometryRepresentation')
bottomDisplay.Representation = 'Surface'
bottomDisplay.ColorArrayName = [None, '']
bottomDisplay.OSPRayScaleArray = 'Normals'
bottomDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
bottomDisplay.OSPRayMaterial = 'scratched'

back = Plane()
back.Origin = [-10.0, -10.0, -10.0]
back.Point1 = [18.0, -10.0, -10.0]
back.Point2 = [-10.0, 10.0, -10.0]
backDisplay = Show(back, renderView, 'GeometryRepresentation')
backDisplay.Representation = 'Surface'
backDisplay.ColorArrayName = [None, '']
backDisplay.AmbientColor = [0.7] * 3
backDisplay.DiffuseColor = [0.7] * 3


anim = GetAnimationScene()
timekeeper = anim.TimeKeeper
pattern = "a_{:}.png"
basenames = list(map(os.path.basename, args.files))
steps = [re.findall("_([0-9]*)", f)[0] for f in basenames]
for i,time in enumerate(timekeeper.TimestepValues):
    fname = pattern.format(steps[i])
    if os.path.isfile(fname) and not args.force:
        Log("skip existing {:}".format(fname))
        continue
    timekeeper.Time = time
    Log("{:}/{:}: {:}".format(i + 1, len(steps), fname))
    SaveScreenshot(fname, renderView)

