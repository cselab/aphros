# Python script to generate images, using ParaView/5.6.0-CrayGNU-18.08-OSMesa
# the command used is SaveAnimation(imgfilename, FrameWindow=[start, end])
# and it requires threee argument given about the pvbatch script.py command
#
# Tesed Tue Dec  4 09:21:55 CET 2018

import sys
import os
if (len(sys.argv) != 4):
  print("missing arguments: pgname imgfile start-index end-index")
  exit()

imgfile = str(sys.argv[1])
start = int(sys.argv[2])
end = int(sys.argv[3])
#W = 6000
W = 1000
q = 2**0.5
H = int(W*q)
SamplesPerPixel = 200


print("rendering images ", imgfile, " from index ", start, " to index ", end, " at ", SamplesPerPixel, " samples per pixel")

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Light'
light9 = CreateLight()
light9.Intensity = 50.0
light9.Type = 'Positional'
light9.Position = [0.5, 10.0, 0.5]
light9.FocalPoint = [0.5, 2.0, 0.5]
light9.Radius = 2.0

# Create a new 'Light'
light10 = CreateLight()
light10.Intensity = 200.0
light10.Type = 'Positional'
light10.Position = [10.0, 25.0, 10.0]
light10.FocalPoint = [0.5, 2.0, 0.5]
light10.ConeAngle = 180.0
light10.Radius = 2.0

mf = "m.json"
mf = os.path.abspath(os.path.realpath(mf))

if not os.path.isfile(mf):
    print("writing {:}".format(mf))
    open(mf, 'w').write('''
    {
      "family" : "OSPRay",
      "version" : "0.0",
      "materials" : {
        "water5" : {
          "type": "Glass",
          "doubles" : {
              "attenuationColor" : [0.22, 0.34, 0.47],
              "attenuationDistance" : [5.],
              "etaInside" : [1.33]
          }
        }
      }
    }
    ''')
          #"attenuationColor" : [0.22, 0.34, 0.47],
          #"attenuationColor" : [0.2, 0.5, 0.8],
          #"attenuationDistance" : [3.],
          #"etaInside" : [1.2]

materialLibrary1 = GetMaterialLibrary()
print("using materials: {:}".format(mf))
materialLibrary1.LoadMaterials = mf

# create light
# create light
# create light
# Create a new 'Render View'
renderView1 = GetRenderView()
#renderView1.ViewSize = [1920, 1080]



renderView1.ViewSize = [W, H]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5, 2.0, 0.5]
renderView1.UseLight = 0
renderView1.StereoType = 0
#renderView1.CameraPosition = [3.28496, 4.4988399999999995, 8.583929999999999]
#renderView1.CameraFocalPoint = [0.7225109999999996, 2.199650000000001, 1.14589]
#renderView1.CameraViewUp = [-0.10387299335602392, 0.9597639386111011, -0.260889983312825]
renderView1.CameraPosition = [1.973490173381768, 5.580715746106084, 8.31740378951455]
renderView1.CameraFocalPoint = [0.5982548420596657, 2.2742135647652844, 0.9449801931947877]
renderView1.CameraViewUp = [-0.07674790436792649, 0.9150106529035253, -0.39606219744908255]
renderView1.CameraParallelScale = 2.12132
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.EnableOSPRay = 1
renderView1.OSPRayRenderer = 'pathtracer'
renderView1.AmbientSamples = 1
renderView1.SamplesPerPixel = SamplesPerPixel
renderView1.AdditionalLights = [light9, light10]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Box'
box1 = Box()
box1.XLength = 100.0
box1.YLength = 0.05
box1.ZLength = 100.0
box1.Center = [0.5, -0.024, 0.5]

# create a list of pathnames to feed to the reader
import glob
fnames = glob.glob("vf*.xmf")
fnames.sort()

reader = XDMFReader(FileNames=fnames)
reader.CellArrayStatus = ['vf']
reader.GridStatus = ['Grid_2']
reader.UpdatePipeline()
vf_xmf = reader

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=vf_xmf)

# create a new 'Iso Volume'
isoVolume1 = IsoVolume(Input=cellDatatoPointData1)
isoVolume1.InputScalars = ['POINTS', 'vf']
isoVolume1.ThresholdRange = [0.0, 0.5]

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(Input=isoVolume1)

# create a new 'Generate Surface Normals'
generateSurfaceNormals1 = GenerateSurfaceNormals(Input=extractSurface1)
generateSurfaceNormals1.FeatureAngle = 80.0

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 7
cylinder1.Height = 11.0
cylinder1.Radius = 50.0
cylinder1.Center = [0.5, 5.0, 0.5]
cylinder1.Capping = 0

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from vf_xmf
vf_xmfDisplay = Show(vf_xmf, renderView1)

# trace defaults for the display properties.
vf_xmfDisplay.Representation = 'Outline'
vf_xmfDisplay.AmbientColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.ColorArrayName = ['CELLS', '']
vf_xmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
vf_xmfDisplay.SelectOrientationVectors = 'None'
vf_xmfDisplay.ScaleFactor = 0.4
vf_xmfDisplay.SelectScaleArray = 'vf'
vf_xmfDisplay.GlyphType = 'Arrow'
vf_xmfDisplay.GlyphTableIndexArray = 'vf'
vf_xmfDisplay.GaussianRadius = 0.02
vf_xmfDisplay.SetScaleArray = [None, '']
vf_xmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
vf_xmfDisplay.OpacityArray = [None, '']
vf_xmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
vf_xmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
vf_xmfDisplay.SelectionCellLabelFontFile = ''
vf_xmfDisplay.SelectionPointLabelFontFile = ''
vf_xmfDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
vf_xmfDisplay.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
vf_xmfDisplay.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
vf_xmfDisplay.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
vf_xmfDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.DataAxesGrid.XTitleFontFile = ''
vf_xmfDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.DataAxesGrid.YTitleFontFile = ''
vf_xmfDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.DataAxesGrid.ZTitleFontFile = ''
vf_xmfDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.DataAxesGrid.XLabelFontFile = ''
vf_xmfDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.DataAxesGrid.YLabelFontFile = ''
vf_xmfDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
vf_xmfDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.PolarAxes.PolarAxisTitleFontFile = ''
vf_xmfDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.PolarAxes.PolarAxisLabelFontFile = ''
vf_xmfDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
vf_xmfDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
vf_xmfDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from cellDatatoPointData1
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Outline'
cellDatatoPointData1Display.AmbientColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.ColorArrayName = ['POINTS', '']
cellDatatoPointData1Display.OSPRayScaleArray = 'vf'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 0.4
cellDatatoPointData1Display.SelectScaleArray = 'vf'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'vf'
cellDatatoPointData1Display.GaussianRadius = 0.02
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'vf']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'vf']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.SelectionCellLabelFontFile = ''
cellDatatoPointData1Display.SelectionPointLabelFontFile = ''
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cellDatatoPointData1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cellDatatoPointData1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cellDatatoPointData1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
cellDatatoPointData1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.XTitleFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.YTitleFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.ZTitleFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.XLabelFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.YLabelFontFile = ''
cellDatatoPointData1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
cellDatatoPointData1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.PolarAxes.PolarAxisTitleFontFile = ''
cellDatatoPointData1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.PolarAxes.PolarAxisLabelFontFile = ''
cellDatatoPointData1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.PolarAxes.LastRadialAxisTextFontFile = ''
cellDatatoPointData1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from box1
box1Display = Show(box1, renderView1)

# trace defaults for the display properties.
box1Display.Representation = 'Surface'
box1Display.AmbientColor = [0.0, 0.0, 0.0]
box1Display.ColorArrayName = ['POINTS', '']
box1Display.OSPRayScaleArray = 'Normals'
box1Display.OSPRayScaleFunction = 'PiecewiseFunction'
box1Display.SelectOrientationVectors = 'None'
box1Display.ScaleFactor = 0.1
box1Display.SelectScaleArray = 'None'
box1Display.GlyphType = 'Arrow'
box1Display.GlyphTableIndexArray = 'None'
box1Display.GaussianRadius = 0.005
box1Display.SetScaleArray = ['POINTS', 'Normals']
box1Display.ScaleTransferFunction = 'PiecewiseFunction'
box1Display.OpacityArray = ['POINTS', 'Normals']
box1Display.OpacityTransferFunction = 'PiecewiseFunction'
box1Display.DataAxesGrid = 'GridAxesRepresentation'
box1Display.SelectionCellLabelFontFile = ''
box1Display.SelectionPointLabelFontFile = ''
box1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
box1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.0, 0.9852216839790344, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
box1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
box1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
box1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
box1Display.DataAxesGrid.XTitleFontFile = ''
box1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
box1Display.DataAxesGrid.YTitleFontFile = ''
box1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
box1Display.DataAxesGrid.ZTitleFontFile = ''
box1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
box1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
box1Display.DataAxesGrid.XLabelFontFile = ''
box1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
box1Display.DataAxesGrid.YLabelFontFile = ''
box1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
box1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
box1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
box1Display.PolarAxes.PolarAxisTitleFontFile = ''
box1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
box1Display.PolarAxes.PolarAxisLabelFontFile = ''
box1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
box1Display.PolarAxes.LastRadialAxisTextFontFile = ''
box1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
box1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from isoVolume1
isoVolume1Display = Show(isoVolume1, renderView1)

# get color transfer function/color map for 'vf'
vfLUT = GetColorTransferFunction('vf')
vfLUT.RGBPoints = [0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0]
vfLUT.ColorSpace = 'RGB'
vfLUT.NanColor = [1.0, 0.0, 0.0]
vfLUT.Discretize = 0
vfLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'vf'
vfPWF = GetOpacityTransferFunction('vf')
vfPWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
vfPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
isoVolume1Display.Representation = 'Surface'
isoVolume1Display.AmbientColor = [0.0, 0.0, 0.0]
isoVolume1Display.ColorArrayName = ['CELLS', 'vf']
isoVolume1Display.LookupTable = vfLUT
isoVolume1Display.OSPRayScaleFunction = 'PiecewiseFunction'
isoVolume1Display.OSPRayMaterial = 'glass'
isoVolume1Display.SelectOrientationVectors = 'None'
isoVolume1Display.ScaleFactor = 0.31875000000000003
isoVolume1Display.SelectScaleArray = 'vf'
isoVolume1Display.GlyphType = 'Arrow'
isoVolume1Display.GlyphTableIndexArray = 'vf'
isoVolume1Display.GaussianRadius = 0.0159375
isoVolume1Display.SetScaleArray = [None, '']
isoVolume1Display.ScaleTransferFunction = 'PiecewiseFunction'
isoVolume1Display.OpacityArray = [None, '']
isoVolume1Display.OpacityTransferFunction = 'PiecewiseFunction'
isoVolume1Display.DataAxesGrid = 'GridAxesRepresentation'
isoVolume1Display.SelectionCellLabelFontFile = ''
isoVolume1Display.SelectionPointLabelFontFile = ''
isoVolume1Display.PolarAxes = 'PolarAxesRepresentation'
isoVolume1Display.ScalarOpacityFunction = vfPWF
isoVolume1Display.ScalarOpacityUnitDistance = 0.3361874622178882

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
isoVolume1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
isoVolume1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
isoVolume1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
isoVolume1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
isoVolume1Display.DataAxesGrid.XTitleFontFile = ''
isoVolume1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
isoVolume1Display.DataAxesGrid.YTitleFontFile = ''
isoVolume1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
isoVolume1Display.DataAxesGrid.ZTitleFontFile = ''
isoVolume1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
isoVolume1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
isoVolume1Display.DataAxesGrid.XLabelFontFile = ''
isoVolume1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
isoVolume1Display.DataAxesGrid.YLabelFontFile = ''
isoVolume1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
isoVolume1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
isoVolume1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
isoVolume1Display.PolarAxes.PolarAxisTitleFontFile = ''
isoVolume1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
isoVolume1Display.PolarAxes.PolarAxisLabelFontFile = ''
isoVolume1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
isoVolume1Display.PolarAxes.LastRadialAxisTextFontFile = ''
isoVolume1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
isoVolume1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from extractSurface1
extractSurface1Display = Show(extractSurface1, renderView1)

# trace defaults for the display properties.
extractSurface1Display.Representation = 'Surface'
extractSurface1Display.AmbientColor = [0.0, 0.0, 0.0]
extractSurface1Display.ColorArrayName = ['CELLS', 'vf']
extractSurface1Display.LookupTable = vfLUT
extractSurface1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractSurface1Display.OSPRayMaterial = 'glass'
extractSurface1Display.SelectOrientationVectors = 'None'
extractSurface1Display.ScaleFactor = 0.4
extractSurface1Display.SelectScaleArray = 'vf'
extractSurface1Display.GlyphType = 'Arrow'
extractSurface1Display.GlyphTableIndexArray = 'vf'
extractSurface1Display.GaussianRadius = 0.02
extractSurface1Display.SetScaleArray = [None, '']
extractSurface1Display.ScaleTransferFunction = 'PiecewiseFunction'
extractSurface1Display.OpacityArray = [None, '']
extractSurface1Display.OpacityTransferFunction = 'PiecewiseFunction'
extractSurface1Display.DataAxesGrid = 'GridAxesRepresentation'
extractSurface1Display.SelectionCellLabelFontFile = ''
extractSurface1Display.SelectionPointLabelFontFile = ''
extractSurface1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
extractSurface1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extractSurface1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extractSurface1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
extractSurface1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
extractSurface1Display.DataAxesGrid.XTitleFontFile = ''
extractSurface1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
extractSurface1Display.DataAxesGrid.YTitleFontFile = ''
extractSurface1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
extractSurface1Display.DataAxesGrid.ZTitleFontFile = ''
extractSurface1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
extractSurface1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
extractSurface1Display.DataAxesGrid.XLabelFontFile = ''
extractSurface1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
extractSurface1Display.DataAxesGrid.YLabelFontFile = ''
extractSurface1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
extractSurface1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
extractSurface1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
extractSurface1Display.PolarAxes.PolarAxisTitleFontFile = ''
extractSurface1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
extractSurface1Display.PolarAxes.PolarAxisLabelFontFile = ''
extractSurface1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
extractSurface1Display.PolarAxes.LastRadialAxisTextFontFile = ''
extractSurface1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
extractSurface1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from generateSurfaceNormals1
generateSurfaceNormals1Display = Show(generateSurfaceNormals1, renderView1)

# trace defaults for the display properties.
generateSurfaceNormals1Display.Representation = 'Surface'
generateSurfaceNormals1Display.AmbientColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.ColorArrayName = ['POINTS', '']
generateSurfaceNormals1Display.OSPRayScaleArray = 'Normals'
generateSurfaceNormals1Display.OSPRayScaleFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.OSPRayMaterial = 'water5'
generateSurfaceNormals1Display.SelectOrientationVectors = 'None'
generateSurfaceNormals1Display.ScaleFactor = 0.4
generateSurfaceNormals1Display.SelectScaleArray = 'vf'
generateSurfaceNormals1Display.GlyphType = 'Arrow'
generateSurfaceNormals1Display.GlyphTableIndexArray = 'vf'
generateSurfaceNormals1Display.GaussianRadius = 0.02
generateSurfaceNormals1Display.SetScaleArray = ['POINTS', 'Normals']
generateSurfaceNormals1Display.ScaleTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.OpacityArray = ['POINTS', 'Normals']
generateSurfaceNormals1Display.OpacityTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.DataAxesGrid = 'GridAxesRepresentation'
generateSurfaceNormals1Display.SelectionCellLabelFontFile = ''
generateSurfaceNormals1Display.SelectionPointLabelFontFile = ''
generateSurfaceNormals1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
generateSurfaceNormals1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
generateSurfaceNormals1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
generateSurfaceNormals1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
generateSurfaceNormals1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.DataAxesGrid.XTitleFontFile = ''
generateSurfaceNormals1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.DataAxesGrid.YTitleFontFile = ''
generateSurfaceNormals1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.DataAxesGrid.ZTitleFontFile = ''
generateSurfaceNormals1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.DataAxesGrid.XLabelFontFile = ''
generateSurfaceNormals1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.DataAxesGrid.YLabelFontFile = ''
generateSurfaceNormals1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
generateSurfaceNormals1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.PolarAxes.PolarAxisTitleFontFile = ''
generateSurfaceNormals1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.PolarAxes.PolarAxisLabelFontFile = ''
generateSurfaceNormals1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.PolarAxes.LastRadialAxisTextFontFile = ''
generateSurfaceNormals1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
generateSurfaceNormals1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from cylinder1
cylinder1Display = Show(cylinder1, renderView1)

# trace defaults for the display properties.
cylinder1Display.Representation = 'Surface'
cylinder1Display.AmbientColor = [0.0, 0.0, 0.0]
cylinder1Display.ColorArrayName = [None, '']
cylinder1Display.OSPRayScaleArray = 'Normals'
cylinder1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cylinder1Display.SelectOrientationVectors = 'None'
cylinder1Display.ScaleFactor = 0.1
cylinder1Display.SelectScaleArray = 'None'
cylinder1Display.GlyphType = 'Arrow'
cylinder1Display.GlyphTableIndexArray = 'None'
cylinder1Display.GaussianRadius = 0.005
cylinder1Display.SetScaleArray = ['POINTS', 'Normals']
cylinder1Display.ScaleTransferFunction = 'PiecewiseFunction'
cylinder1Display.OpacityArray = ['POINTS', 'Normals']
cylinder1Display.OpacityTransferFunction = 'PiecewiseFunction'
cylinder1Display.DataAxesGrid = 'GridAxesRepresentation'
cylinder1Display.SelectionCellLabelFontFile = ''
cylinder1Display.SelectionPointLabelFontFile = ''
cylinder1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cylinder1Display.OSPRayScaleFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cylinder1Display.ScaleTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cylinder1Display.OpacityTransferFunction.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
cylinder1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
cylinder1Display.DataAxesGrid.XTitleFontFile = ''
cylinder1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
cylinder1Display.DataAxesGrid.YTitleFontFile = ''
cylinder1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
cylinder1Display.DataAxesGrid.ZTitleFontFile = ''
cylinder1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
cylinder1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
cylinder1Display.DataAxesGrid.XLabelFontFile = ''
cylinder1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
cylinder1Display.DataAxesGrid.YLabelFontFile = ''
cylinder1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
cylinder1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
cylinder1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
cylinder1Display.PolarAxes.PolarAxisTitleFontFile = ''
cylinder1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
cylinder1Display.PolarAxes.PolarAxisLabelFontFile = ''
cylinder1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
cylinder1Display.PolarAxes.LastRadialAxisTextFontFile = ''
cylinder1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
cylinder1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# setup the color legend parameters for each legend in this view

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.InterpretValuesAsCategories = 1
vtkBlockColorsLUT.AnnotationsInitialized = 1
vtkBlockColorsLUT.RGBPoints = [0.005, 0.231373, 0.298039, 0.752941, 0.0125, 0.865003, 0.865003, 0.865003, 0.02, 0.705882, 0.0156863, 0.14902]
vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
vtkBlockColorsLUT.ActiveAnnotatedValues = ['0', '1']
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.6299992370489051, 0.6299992370489051, 1.0, 0.6699931334401464, 0.5000076295109483, 0.3300068665598535, 1.0, 0.5000076295109483, 0.7499961852445258, 0.5300068665598535, 0.3499961852445258, 0.7000076295109483, 1.0, 0.7499961852445258, 0.5000076295109483]
vtkBlockColorsLUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# get color legend/bar for vtkBlockColorsLUT in view renderView1
vtkBlockColorsLUTColorBar = GetScalarBar(vtkBlockColorsLUT, renderView1)
vtkBlockColorsLUTColorBar.Title = 'vtkBlockColors'
vtkBlockColorsLUTColorBar.ComponentTitle = ''
vtkBlockColorsLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
vtkBlockColorsLUTColorBar.TitleFontFile = ''
vtkBlockColorsLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
vtkBlockColorsLUTColorBar.LabelFontFile = ''

# set color bar visibility
vtkBlockColorsLUTColorBar.Visibility = 0

# get color legend/bar for vfLUT in view renderView1
vfLUTColorBar = GetScalarBar(vfLUT, renderView1)
vfLUTColorBar.Title = 'vf'
vfLUTColorBar.ComponentTitle = ''
vfLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
vfLUTColorBar.TitleFontFile = ''
vfLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
vfLUTColorBar.LabelFontFile = ''

# set color bar visibility
vfLUTColorBar.Visibility = 0

# hide data in view
Hide(vf_xmf, renderView1)

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# hide data in view
Hide(isoVolume1, renderView1)

# hide data in view
Hide(extractSurface1, renderView1)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')
vtkBlockColorsPWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------

SaveScreenshot(imgfile)
