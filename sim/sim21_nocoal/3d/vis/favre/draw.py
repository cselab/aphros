#!/usr/bin/env pvbatch

# state file generated using paraview version 5.6.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000, 1000]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [2.864659435788656, 1.400169527753375, -1.1935958157385422]
renderView1.CameraFocalPoint = [0.7618904044343615, 0.5023327450596681, 0.36183797116664934]
renderView1.CameraViewUp = [-0.26122254713626997, 0.9458257135082125, 0.19281208606910466]
renderView1.CameraParallelScale = 0.8660254037844386
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableOSPRay = 1
renderView1.Shadows = 1
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 10
renderView1.ProgressivePasses = 1
renderView1.OSPRayMaterialLibrary = materialLibrary1


s_0500vtk = LegacyVTKReader(FileNames=['s_0500.vtk'])

vf_0500xmf = XDMFReader(FileNames=['vf_0500.xmf'])
vf_0500xmf.CellArrayStatus = ['vf']
vf_0500xmf.GridStatus = ['Grid_2']

cellDatatoPointData1 = CellDatatoPointData(Input=vf_0500xmf)

contour1 = Contour(Input=cellDatatoPointData1)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'


cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)
cellDatatoPointData1Display.Representation = 'Outline'
cellDatatoPointData1Display.AmbientColor = [0.0, 0.0, 0.0]
cellDatatoPointData1Display.LineWidth = 1

contour1Display = Show(contour1, renderView1)
contour1Display.Representation = 'Surface'
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.Opacity = 0.5

s_0500vtkDisplay = Show(s_0500vtk, renderView1)
s_0500vtkDisplay.Representation = 'Surface'
s_0500vtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
s_0500vtkDisplay.Opacity = 0.5


fn = "polygons.png"
Hide(s_0500vtk, renderView1)
Show(contour1, renderView1)
SaveScreenshot(fn, renderView1)
print(fn)

fn = "contour.png"
Show(s_0500vtk, renderView1)
Hide(contour1, renderView1)
SaveScreenshot(fn, renderView1)
print(fn)

