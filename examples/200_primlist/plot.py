#!/usr/bin/env pvbatch

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1080, 1080]
renderView1.OrientationAxesVisibility = 1
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CameraPosition = [2.850257956930087, 1.3875069244683658, 1.2981201068256059]
renderView1.CameraFocalPoint = [-0.062126483577120056, 0.21088117807884124, 0.14498475209410736]
renderView1.CameraViewUp = [-0.31233559104266445, -0.14673926220996603, 0.9385702251265509]
renderView1.Background = [1.0, 1.0, 1.0]

sm = LegacyVTKReader(FileNames=['sm_0000.vtk'])

surf = Threshold(Input=sm)
surf.Scalars = ['CELLS', 'cl']
surf.ThresholdRange = [0.0, 10]

surfDisplay = Show(surf, renderView1)
clLUT = GetColorTransferFunction('cl')
clLUT.AutomaticRescaleRangeMode = 'Never'
clLUT.RGBPoints = [
        -1.5, 0.0, 0.0, 0.5625,
        -0.5, 0.0, 0.0, 1.0,
        0.3, 0.0, 1.0, 1.0,
        1.8, 0.5, 1.0, 0.5,
        3.5, 1.0, 1.0, 0.0,
        5, 1, 0.0, 0.0,
        ]
clLUT.ColorSpace = 'RGB'
clLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
clLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
surfDisplay.Representation = 'Surface'
surfDisplay.ColorArrayName = ['CELLS', 'cl']
surfDisplay.LookupTable = clLUT

fn = "a.jpg"
print(fn)
SaveScreenshot(fn)
