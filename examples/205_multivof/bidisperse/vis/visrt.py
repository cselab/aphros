#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np
import paratools

parser = argparse.ArgumentParser(
    description="Pathtracer rendering of bubbles and boundaries.")
parser.add_argument('files', nargs='+', help="list of data files 'sm_*.vtk'")
parser.add_argument('--force',
                    action="store_true",
                    help="overwrite existing files")
parser.add_argument('--cover',
                    action="store_true",
                    help="use settings for cover")
parser.add_argument('--draft', action="store_true", help="less samples")
parser.add_argument('--resolution',
                    nargs=2,
                    type=int,
                    default=[1920, 1080],
                    help="image resolution")
args = parser.parse_args()

source_bcvtk = LegacyVTKReader(
    FileNames=paratools.ReplaceFilename(args.files[:1], "bc.vtk"))
source_sm = LegacyVTKReader(FileNames=args.files)
sources_ft, timearrays = paratools.ApplyForceTime([source_sm])
source_sm, = sources_ft

light1 = CreateLight()
light1.Intensity = 8.0
light1.Type = 'Positional'
light1.Position = [
    0.006431806423015062, 1.4778652489672959, 1.6991249052313981
]
light1.FocalPoint = [
    1.9781192660676643, -0.1963031716820185, -1.997106578269941
]
#light1.DiffuseColor = [0.95, 0.75, 0.95]
light1.DiffuseColor = [1, 1, 1]
light1.ConeAngle = 27.601142162393426
light1.Radius = 2.0

light2 = CreateLight()
light2.Intensity = 8.0
light2.Type = 'Positional'
light2.Position = [2.7913862587239078, 0.21806733382668228, 1.3226099304191168]
light2.FocalPoint = [
    -0.08931695523427821, 0.6885663380128587, -2.1172589913563358
]
#light2.DiffuseColor = [0.8, 1.0, 0.8]
light2.DiffuseColor = [1, 1, 1]
light2.Radius = 2.0

light3 = CreateLight()
light3.Coords = 'Ambient'
light3.Intensity = 0.5

light4 = CreateLight()
light4.Intensity = 10.0
light4.Type = 'Positional'
light4.Position = [0.5153479625151784, -1.050258016384658, 1.4989357787750535]
light4.FocalPoint = [1.2467274431718631, 2.1, -1.529410314885414]
light4.ConeAngle = 10.0
light4.Radius = 2.0

mf = "m.json"
open(mf, 'w').write('''
{
  "family" : "OSPRay",
  "version" : "0.0",
  "materials" : {
    "waterwhite" : {
      "type": "Glass",
      "doubles" : {
          "attenuationColor" : [0.35, 0.35, 0.35],
          "attenuationDistance" : [0.5],
          "eta" : [1.33]
      }
    }
  }
}
''')
materialLibrary1 = GetMaterialLibrary()
materialLibrary1.LoadMaterials = mf

renderView1 = CreateView('RenderView')
renderView1.ViewSize = args.resolution
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
renderView1.KeyLightWarmth = 0.5
renderView1.FillLightWarmth = 0.5
if args.cover:
    renderView1.CameraPosition = [0.58541, 0.055019, 0.58551]
    renderView1.CameraFocalPoint = [2.8025, 5.2331, -4.2809]
    renderView1.CameraViewUp = [0.35988, 0.55249, 0.75183]
else:
    renderView1.CameraPosition = [
        0.5450666131665406, 0.26691967145900836, 0.3104520152186882
    ]
    renderView1.CameraFocalPoint = [
        2.7119908810557116, 2.7184911650856627, -2.7954705446201813
    ]
    renderView1.CameraViewUp = [
        0.5120912672272875, 0.4635107901442511, 0.7231322710606981
    ]
renderView1.CameraViewAngle = 60.0
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.1676301748696558
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay pathtracer'
renderView1.SamplesPerPixel = 5 if args.draft else 100
renderView1.AdditionalLights = [light1, light2, light3, light4]
renderView1.OSPRayMaterialLibrary = materialLibrary1

ebvf = XDMFReader(FileNames=['../ebvf_0000.xmf', '../ebvf_0001.xmf'])
ebvf.CellArrayStatus = ['ebvf']
ebvf.GridStatus = ['Grid_0']
celltopoint = CellDatatoPointData(Input=ebvf)
contouter = Contour(Input=celltopoint)
contouter.ContourBy = ['POINTS', 'ebvf']
contouter.Isosurfaces = [0.25]

continner = Contour(Input=celltopoint)
continner.ContourBy = ['POINTS', 'ebvf']
continner.Isosurfaces = [0.3]

clip = Clip(Input=contouter)
clip.ClipType = 'Plane'
clip.HyperTreeGridClipper = 'Plane'
clip.Scalars = ['POINTS', 'ebvf']
clip.ClipType.Origin = [1.0049999952316284, 0.5150625295937061, 0.044]
clip.ClipType.Normal = [0.0, 0.0, 1.0]

bubnn = Calculator(Input=source_sm)
bubnn.ResultNormals = 1
bubnn.ResultArrayName = 'Normals'
bubnn.Function = 'nn'

appendgeom = AppendGeometry(Input=[bubnn, continner])
gensurfnorm = GenerateSurfaceNormals(Input=appendgeom)

box = paratools.GetBoundingBox(contouter)

bottom = Plane()
bottom_z = box[0][2] - 0.001
bottom.Origin = [-3.0, -3.0, bottom_z]
bottom.Point1 = [6.0, -3.0, bottom_z]
bottom.Point2 = [-3.0, 6.0, bottom_z]
bottomDisplay = Show(bottom, renderView1, 'GeometryRepresentation')
bottomDisplay.Representation = 'Surface'
bottomDisplay.AmbientColor = [0.85] * 3
bottomDisplay.ColorArrayName = [None, '']
bottomDisplay.DiffuseColor = [0.85] * 3

gensurfnormDisplay = Show(gensurfnorm, renderView1, 'GeometryRepresentation')
gensurfnormDisplay.Representation = 'Surface'
gensurfnormDisplay.ColorArrayName = [None, '']
gensurfnormDisplay.OSPRayMaterial = 'waterwhite'

clipDisplay = Show(clip, renderView1, 'UnstructuredGridRepresentation')
clipDisplay.Representation = 'Surface'
clipDisplay.ColorArrayName = ['POINTS', '']
clipDisplay.Representation = 'Surface'
clipDisplay.AmbientColor = [0.85] * 3
clipDisplay.ColorArrayName = [None, '']
clipDisplay.DiffuseColor = [0.85] * 3

steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)
