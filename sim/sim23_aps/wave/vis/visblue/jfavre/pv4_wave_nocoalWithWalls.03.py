Batch = True

if Batch:
  import sys

  if (len(sys.argv) != 4):
    print("missing arguments: pgname imgfile start-index end-index")
    exit()

  imgfile = str(sys.argv[1])
  start = int(sys.argv[2])
  end = int(sys.argv[3])

if Batch:
  ProgressivePasses = 0
  SamplesPerPixel = 200
  ViewSize = [1920, 1080]
else:
  ViewSize = [1280, 720]
  ProgressivePasses = 10000
  SamplesPerPixel = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

mf1 = "/users/jfavre/Projects/ParaView/ospray_mats.json"
materialLibrary1 = GetMaterialLibrary()
print("using materials: {:}".format(mf1))
materialLibrary1.LoadMaterials = mf1

light2 = CreateLight()
light2.Intensity = 60.0
light2.Type = 'Positional'
light2.Position = [2.0, 8.0, 1.0]
light2.FocalPoint = [1.0, 1.0, 0.5]
light2.ConeAngle = 11.0
light2.Radius = 0.5

light11 = CreateLight()
light11.Intensity = 100.0
light11.Type = 'Positional'
light11.Position = [2.0, 20.0, 2.0]
light11.FocalPoint = [1.0, 0.0, 0.5]
light11.ConeAngle = 23.0
light11.Radius = 5.0


# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1666, 968]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CenterOfRotation = [0.0, 0.3267297148704529, 0.0]
renderView1.UseLight = 0

#renderView1.CameraPosition = [2.1509697331816295, 1.784817137207785, 3.199136175464221]
#renderView1.CameraFocalPoint = [-17.013224735263993, -23.99284889590715, -41.02167171453066]
#renderView1.CameraViewUp = [-0.26277401995216415, 0.8787720902285784, -0.3983835185766037]
renderView1.CameraPosition = [1.673947199541654, 1.155534385955274, 3.017460244409912]
renderView1.CameraFocalPoint = [0.49676167169496643, -0.25002651930486053, -1.3449860440897603]
renderView1.CameraViewUp = [-0.0893981734829278, 0.9547941520217919, -0.2835067791833273]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 14.14590938422056
renderView1.Background = [0.1803921568627451, 0.20392156862745098, 0.21176470588235294]
renderView1.EnableRayTracing = 1
renderView1.LightScale=1
renderView1.BackEnd = 'OSPRay pathtracer'
renderView1.AmbientSamples = 0
renderView1.SamplesPerPixel = SamplesPerPixel
renderView1.ProgressivePasses = ProgressivePasses
renderView1.AdditionalLights = [light2, light11]
renderView1.OSPRayMaterialLibrary = materialLibrary1
renderView1.OrientationAxesVisibility=0
SetActiveView(None)

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

SetActiveView(renderView1)

if Batch:
  fnames=[format("/scratch/snx3000/karnakov/pub/aps/4_wave_nocoal/nx256l2wang/sm_%04d.vtk" % i) for i in range(0,656, 1)]
else:
  fnames=[format("/mnt/data/Karnakov/aps/4_wave_nocoal/nx256l2wang65/sm_%04d.vtk" % i) for i in range(0,656, 1)]

# create a new 'XDMF Reader'
reader = LegacyVTKReader(FileNames=fnames)
reader.UpdatePipelineInformation()

surf = Calculator(Input=reader)
surf.ResultNormals = 1
surf.AttributeType = 'Point Data'
surf.ResultArrayName = 'Normals'
surf.Function = 'nn'

# create a new 'Plane'
plane1 = Plane()
plane1.Origin = [-20.0, -0.01, -35.0]
plane1.Point1 = [-20.0, -0.01, 5.0]
plane1.Point2 = [20.0, -0.01, -35.0]

plane1Display = Show(plane1, renderView1)
plane1Display.Representation = 'Surface'
plane1Display.OSPRayMaterial = 'checker grey'
plane1Display.ColorArrayName = ['POINTS', '']

bottom = Plane(guiName="bottom")
bottom.Origin = [0.0, 0.00015, 0.0]
bottom.Point1 = [2.0, 0.00015, 0.0]
bottom.Point2 = [0.0, 0.00015, 1.0]
bottomDisplay = Show(bottom, renderView1)
bottomDisplay.Representation = 'Surface'
bottomDisplay.ColorArrayName = [None, '']

# show data from surf
surfDisplay = Show(surf, renderView1)

# trace defaults for the display properties.
surfDisplay.Representation = 'Surface'
surfDisplay.ColorArrayName = ['POINTS', '']
surfDisplay.OSPRayScaleArray = 'Normals'
surfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
surfDisplay.OSPRayMaterial = 'water'
surfDisplay.ScaleFactor = 0.2
surfDisplay.SelectScaleArray = 'c'
surfDisplay.GlyphType = 'Arrow'
surfDisplay.GlyphTableIndexArray = 'c'
surfDisplay.GaussianRadius = 0.01
surfDisplay.SetScaleArray = ['POINTS', 'Normals']
surfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
surfDisplay.OpacityArray = ['POINTS', 'Normals']
surfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
surfDisplay.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
surfDisplay.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.NumberOfFrames = len(reader.TimestepValues)
scene.PlayMode = 'Snap To TimeSteps'

Render()
if Batch:
  SaveAnimation(imgfile, FrameWindow=[start, end], ImageResolution=ViewSize)
