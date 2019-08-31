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
renderView1.ViewSize = [1000,1000]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 0
renderView1.CameraPosition =\
[1.326166076499316, 1.581412512627694, 3.5052432288262967]
renderView1.CameraFocalPoint =\
[0.45722849374168073, 0.4442839115593683, 0.4806707888589165]
renderView1.CameraViewUp =\
[-0.08433283755119116, 0.9404371241304461, -0.3293417496580418]
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.EnableOSPRay = 1
renderView1.Shadows = 1
renderView1.AmbientSamples = 5
renderView1.SamplesPerPixel = 10
renderView1.OSPRayMaterialLibrary = materialLibrary1



def Draw0():
    vfxmf = XDMFReader(FileNames=['vf_0500.xmf'])
    vfxmf.CellArrayStatus = ['vf']
    vfxmf.GridStatus = ['Grid_2']

    cellDatatoPointData1 = CellDatatoPointData(Input=vfxmf)
    contour1 = Contour(Input=cellDatatoPointData1)
    contour1.ContourBy = ['POINTS', 'vf']
    contour1.Isosurfaces = [0.5]
    contour1.PointMergeMethod = 'Uniform Binning'
    contour1Display = Show(contour1, renderView1)
    contour1Display.Representation = 'Surface'
    contour1Display.AmbientColor = [0.0, 0.0, 0.0]

    fn = "contour.png"
    SaveScreenshot(fn, renderView1)
    print(fn)
    Hide(contour1)


def Draw1():
    svtk = LegacyVTKReader(FileNames=['s_0500.vtk'])
    svtkDisplay = Show(svtk, renderView1)
    svtkDisplay.Representation = 'Surface'
    svtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
    svtkDisplay.ColorArrayName = [None, '']

    fn = "polygons.png"
    SaveScreenshot(fn, renderView1)
    print(fn)
    Hide(svtk)


def Draw2():
    smvtk = LegacyVTKReader(FileNames=['sm_0500.vtk'])
    smvtk = GenerateSurfaceNormals(Input=smvtk)
    clip1 = Clip(Input=smvtk)
    clip1.ClipType = 'Plane'
    clip1.Crinkleclip = 1
    clip1.ClipType.Origin = [0.0, 0.99, 0.0]
    clip1.ClipType.Normal = [0.0, 1.0, 0.0]
    surfDisplay = Show(clip1, renderView1)
    surfDisplay.Representation = 'Surface'
    surfDisplay.AmbientColor = [0.0, 0.0, 0.0]
    surfDisplay.ColorArrayName = [None, '']
    fn = "polygons_march.png"
    SaveScreenshot(fn, renderView1)
    print(fn)
    Hide(clip1)

def Draw(case):
    if case == 0: Draw0()
    if case == 1: Draw1()
    if case == 2: Draw2()

import sys
av = sys.argv

if len(av) > 1:
    Draw(int(av[1]))
else:
    for case in [0,1,2]:
        Draw(case)
