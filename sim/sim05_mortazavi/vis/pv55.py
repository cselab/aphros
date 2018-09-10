# state file generated using paraview version 5.5.0-RC3

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.0-RC3

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1548, 897]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [2.0, 0.5, 0.5]
renderView1.StereoType = 0
renderView1.CameraPosition = [-2.1217664365211846, 2.0110134368411257, 7.421328945563331]
renderView1.CameraFocalPoint = [2.0, 0.5, 0.5]
renderView1.CameraViewUp = [0.0918323771390451, 0.9828553088329238, -0.15988200776801023]
renderView1.CameraParallelScale = 1.2671729555822187
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

base = '/home/kpetr/s/euler/s/scratch/sim04_p32a/'

# create a new 'XDMF Reader'
vx_ = XDMFReader(FileNames=[base + '/vx_0.xmf', base + '/vx_1.xmf', base + '/vx_2.xmf', base + '/vx_3.xmf', base + '/vx_4.xmf', base + '/vx_5.xmf', base + '/vx_6.xmf', base + '/vx_7.xmf', base + '/vx_8.xmf', base + '/vx_9.xmf', base + '/vx_10.xmf', base + '/vx_11.xmf', base + '/vx_12.xmf', base + '/vx_13.xmf', base + '/vx_14.xmf', base + '/vx_15.xmf', base + '/vx_16.xmf', base + '/vx_17.xmf', base + '/vx_18.xmf', base + '/vx_19.xmf', base + '/vx_20.xmf', base + '/vx_21.xmf', base + '/vx_22.xmf', base + '/vx_23.xmf', base + '/vx_24.xmf', base + '/vx_25.xmf', base + '/vx_26.xmf', base + '/vx_27.xmf', base + '/vx_28.xmf', base + '/vx_29.xmf', base + '/vx_30.xmf', base + '/vx_31.xmf', base + '/vx_32.xmf', base + '/vx_33.xmf', base + '/vx_34.xmf', base + '/vx_35.xmf', base + '/vx_36.xmf', base + '/vx_37.xmf', base + '/vx_38.xmf', base + '/vx_39.xmf', base + '/vx_40.xmf', base + '/vx_41.xmf', base + '/vx_42.xmf', base + '/vx_43.xmf', base + '/vx_44.xmf', base + '/vx_45.xmf', base + '/vx_46.xmf', base + '/vx_47.xmf', base + '/vx_48.xmf', base + '/vx_49.xmf'])
vx_.CellArrayStatus = ['vx']
vx_.GridStatus = ['Grid_7034']

# create a new 'XDMF Reader'
vy_ = XDMFReader(FileNames=[base + '/vy_0.xmf', base + '/vy_1.xmf', base + '/vy_2.xmf', base + '/vy_3.xmf', base + '/vy_4.xmf', base + '/vy_5.xmf', base + '/vy_6.xmf', base + '/vy_7.xmf', base + '/vy_8.xmf', base + '/vy_9.xmf', base + '/vy_10.xmf', base + '/vy_11.xmf', base + '/vy_12.xmf', base + '/vy_13.xmf', base + '/vy_14.xmf', base + '/vy_15.xmf', base + '/vy_16.xmf', base + '/vy_17.xmf', base + '/vy_18.xmf', base + '/vy_19.xmf', base + '/vy_20.xmf', base + '/vy_21.xmf', base + '/vy_22.xmf', base + '/vy_23.xmf', base + '/vy_24.xmf', base + '/vy_25.xmf', base + '/vy_26.xmf', base + '/vy_27.xmf', base + '/vy_28.xmf', base + '/vy_29.xmf', base + '/vy_30.xmf', base + '/vy_31.xmf', base + '/vy_32.xmf', base + '/vy_33.xmf', base + '/vy_34.xmf', base + '/vy_35.xmf', base + '/vy_36.xmf', base + '/vy_37.xmf', base + '/vy_38.xmf', base + '/vy_39.xmf', base + '/vy_40.xmf', base + '/vy_41.xmf', base + '/vy_42.xmf', base + '/vy_43.xmf', base + '/vy_44.xmf', base + '/vy_45.xmf', base + '/vy_46.xmf', base + '/vy_47.xmf', base + '/vy_48.xmf', base + '/vy_49.xmf'])
vy_.CellArrayStatus = ['vy']
vy_.GridStatus = ['Grid_7037']

# create a new 'XDMF Reader'
vf_ = XDMFReader(FileNames=[base + '/vf_0.xmf', base + '/vf_1.xmf', base + '/vf_2.xmf', base + '/vf_3.xmf', base + '/vf_4.xmf', base + '/vf_5.xmf', base + '/vf_6.xmf', base + '/vf_7.xmf', base + '/vf_8.xmf', base + '/vf_9.xmf', base + '/vf_10.xmf', base + '/vf_11.xmf', base + '/vf_12.xmf', base + '/vf_13.xmf', base + '/vf_14.xmf', base + '/vf_15.xmf', base + '/vf_16.xmf', base + '/vf_17.xmf', base + '/vf_18.xmf', base + '/vf_19.xmf', base + '/vf_20.xmf', base + '/vf_21.xmf', base + '/vf_22.xmf', base + '/vf_23.xmf', base + '/vf_24.xmf', base + '/vf_25.xmf', base + '/vf_26.xmf', base + '/vf_27.xmf', base + '/vf_28.xmf', base + '/vf_29.xmf', base + '/vf_30.xmf', base + '/vf_31.xmf', base + '/vf_32.xmf', base + '/vf_33.xmf', base + '/vf_34.xmf', base + '/vf_35.xmf', base + '/vf_36.xmf', base + '/vf_37.xmf', base + '/vf_38.xmf', base + '/vf_39.xmf', base + '/vf_40.xmf', base + '/vf_41.xmf', base + '/vf_42.xmf', base + '/vf_43.xmf', base + '/vf_44.xmf', base + '/vf_45.xmf', base + '/vf_46.xmf', base + '/vf_47.xmf', base + '/vf_48.xmf', base + '/vf_49.xmf'])
vf_.CellArrayStatus = ['vf']
vf_.GridStatus = ['Grid_7031']

# create a new 'XDMF Reader'
vz_ = XDMFReader(FileNames=[base + '/vz_0.xmf', base + '/vz_1.xmf', base + '/vz_2.xmf', base + '/vz_3.xmf', base + '/vz_4.xmf', base + '/vz_5.xmf', base + '/vz_6.xmf', base + '/vz_7.xmf', base + '/vz_8.xmf', base + '/vz_9.xmf', base + '/vz_10.xmf', base + '/vz_11.xmf', base + '/vz_12.xmf', base + '/vz_13.xmf', base + '/vz_14.xmf', base + '/vz_15.xmf', base + '/vz_16.xmf', base + '/vz_17.xmf', base + '/vz_18.xmf', base + '/vz_19.xmf', base + '/vz_20.xmf', base + '/vz_21.xmf', base + '/vz_22.xmf', base + '/vz_23.xmf', base + '/vz_24.xmf', base + '/vz_25.xmf', base + '/vz_26.xmf', base + '/vz_27.xmf', base + '/vz_28.xmf', base + '/vz_29.xmf', base + '/vz_30.xmf', base + '/vz_31.xmf', base + '/vz_32.xmf', base + '/vz_33.xmf', base + '/vz_34.xmf', base + '/vz_35.xmf', base + '/vz_36.xmf', base + '/vz_37.xmf', base + '/vz_38.xmf', base + '/vz_39.xmf', base + '/vz_40.xmf', base + '/vz_41.xmf', base + '/vz_42.xmf', base + '/vz_43.xmf', base + '/vz_44.xmf', base + '/vz_45.xmf', base + '/vz_46.xmf', base + '/vz_47.xmf', base + '/vz_48.xmf', base + '/vz_49.xmf'])
vz_.CellArrayStatus = ['vz']
vz_.GridStatus = ['Grid_7040']

# create a new 'XDMF Reader'
p_ = XDMFReader(FileNames=[base + '/p_0.xmf', base + '/p_1.xmf', base + '/p_2.xmf', base + '/p_3.xmf', base + '/p_4.xmf', base + '/p_5.xmf', base + '/p_6.xmf', base + '/p_7.xmf', base + '/p_8.xmf', base + '/p_9.xmf', base + '/p_10.xmf', base + '/p_11.xmf', base + '/p_12.xmf', base + '/p_13.xmf', base + '/p_14.xmf', base + '/p_15.xmf', base + '/p_16.xmf', base + '/p_17.xmf', base + '/p_18.xmf', base + '/p_19.xmf', base + '/p_20.xmf', base + '/p_21.xmf', base + '/p_22.xmf', base + '/p_23.xmf',
    base + '/p_24.xmf', base + '/p_25.xmf', base + '/p_26.xmf', base + '/p_27.xmf', base + '/p_28.xmf', base + '/p_29.xmf', base + '/p_30.xmf', base + '/p_31.xmf', base + '/p_32.xmf', base + '/p_33.xmf', base + '/p_34.xmf', base + '/p_35.xmf', base + '/p_36.xmf', base + '/p_37.xmf', base + '/p_38.xmf', base + '/p_39.xmf', base + '/p_40.xmf', base + '/p_41.xmf', base + '/p_42.xmf', base + '/p_43.xmf', base + '/p_44.xmf', base + '/p_45.xmf', base + '/p_46.xmf', base + '/p_47.xmf', base +
    '/p_48.xmf', base + '/p_49.xmf'])
p_.CellArrayStatus = ['p']
p_.GridStatus = ['Grid_7028']

# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(Input=[p_, vf_, vx_, vy_, vz_])

# create a new 'Calculator'
calculator1 = Calculator(Input=appendAttributes1)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'vel'
calculator1.Function = 'vx*iHat+vy*jHat+vz*kHat'

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=calculator1)
resampleToImage1.SamplingDimensions = [256, 64, 64]
resampleToImage1.SamplingBounds = [0.0, 4.0, 0.0, 1.0, 0.0, 1.0]

# create a new 'Gradient Of Unstructured DataSet'
gradientOfUnstructuredDataSet1 = GradientOfUnstructuredDataSet(Input=resampleToImage1)
gradientOfUnstructuredDataSet1.ScalarArray = ['POINTS', 'vel']
gradientOfUnstructuredDataSet1.ComputeGradient = 0
gradientOfUnstructuredDataSet1.ComputeVorticity = 1

# create a new 'Contour'
contour1 = Contour(Input=resampleToImage1)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from calculator1
calculator1Display = Show(calculator1, renderView1)

# trace defaults for the display properties.
calculator1Display.Representation = 'Outline'
calculator1Display.ColorArrayName = ['CELLS', '']
calculator1Display.LineWidth = 3.0
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 0.4
calculator1Display.SelectScaleArray = 'p'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'p'
calculator1Display.GaussianRadius = 0.02
calculator1Display.SetScaleArray = [None, '']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = [None, '']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.SelectionCellLabelFontFile = ''
calculator1Display.SelectionPointLabelFontFile = ''
calculator1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calculator1Display.DataAxesGrid.XTitleFontFile = ''
calculator1Display.DataAxesGrid.YTitleFontFile = ''
calculator1Display.DataAxesGrid.ZTitleFontFile = ''
calculator1Display.DataAxesGrid.XLabelFontFile = ''
calculator1Display.DataAxesGrid.YLabelFontFile = ''
calculator1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calculator1Display.PolarAxes.PolarAxisTitleFontFile = ''
calculator1Display.PolarAxes.PolarAxisLabelFontFile = ''
calculator1Display.PolarAxes.LastRadialAxisTextFontFile = ''
calculator1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from contour1
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']
contour1Display.OSPRayScaleArray = 'Normals'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.30901959836483006
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.GaussianRadius = 0.015450979918241502
contour1Display.SetScaleArray = ['POINTS', 'Normals']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'Normals']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.SelectionCellLabelFontFile = ''
contour1Display.SelectionPointLabelFontFile = ''
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
contour1Display.DataAxesGrid.XTitleFontFile = ''
contour1Display.DataAxesGrid.YTitleFontFile = ''
contour1Display.DataAxesGrid.ZTitleFontFile = ''
contour1Display.DataAxesGrid.XLabelFontFile = ''
contour1Display.DataAxesGrid.YLabelFontFile = ''
contour1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour1Display.PolarAxes.PolarAxisTitleFontFile = ''
contour1Display.PolarAxes.PolarAxisLabelFontFile = ''
contour1Display.PolarAxes.LastRadialAxisTextFontFile = ''
contour1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from gradientOfUnstructuredDataSet1
gradientOfUnstructuredDataSet1Display = Show(gradientOfUnstructuredDataSet1, renderView1)

# get color transfer function/color map for 'Vorticity'
vorticityLUT = GetColorTransferFunction('Vorticity')
vorticityLUT.AutomaticRescaleRangeMode = 'Never'
vorticityLUT.RGBPoints = [0.0, 0.0416667, 0.0, 0.0, 1.2698399999999999, 0.208333, 0.0, 0.0, 2.5396799999999997, 0.375, 0.0, 0.0, 3.8095199999999996, 0.541667, 0.0, 0.0, 5.079370000000001, 0.708333, 0.0, 0.0, 6.349210000000001, 0.854137, 0.0, 0.0, 7.6190500000000005, 0.937488, 0.039062, 0.0, 8.88889, 1.0, 0.208333, 0.0, 10.15873, 1.0, 0.375, 0.0, 11.42857, 1.0, 0.541667, 0.0, 12.698409999999999, 1.0, 0.708333, 0.0, 13.96825, 1.0, 0.858805, 0.03125, 15.238100000000001, 1.0, 0.947392, 0.15625, 16.507939999999998, 1.0, 1.0, 0.3125, 17.77778, 1.0, 1.0, 0.5625, 19.04762, 1.0, 1.0, 0.8125, 20.0, 1.0, 1.0, 1.0]
vorticityLUT.ColorSpace = 'Lab'
vorticityLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Vorticity'
vorticityPWF = GetOpacityTransferFunction('Vorticity')
vorticityPWF.Points = [0.0, 0.0, 0.5, 0.0, 20.0, 0.5197368264198303, 0.5, 0.0]
vorticityPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
gradientOfUnstructuredDataSet1Display.Representation = 'Volume'
gradientOfUnstructuredDataSet1Display.ColorArrayName = ['POINTS', 'Vorticity']
gradientOfUnstructuredDataSet1Display.LookupTable = vorticityLUT
gradientOfUnstructuredDataSet1Display.OSPRayScaleArray = 'Vorticity'
gradientOfUnstructuredDataSet1Display.OSPRayScaleFunction = 'PiecewiseFunction'
gradientOfUnstructuredDataSet1Display.SelectOrientationVectors = 'None'
gradientOfUnstructuredDataSet1Display.ScaleFactor = 0.4
gradientOfUnstructuredDataSet1Display.SelectScaleArray = 'None'
gradientOfUnstructuredDataSet1Display.GlyphType = 'Arrow'
gradientOfUnstructuredDataSet1Display.GlyphTableIndexArray = 'None'
gradientOfUnstructuredDataSet1Display.GaussianRadius = 0.02
gradientOfUnstructuredDataSet1Display.SetScaleArray = ['POINTS', 'Vorticity']
gradientOfUnstructuredDataSet1Display.ScaleTransferFunction = 'PiecewiseFunction'
gradientOfUnstructuredDataSet1Display.OpacityArray = ['POINTS', 'Vorticity']
gradientOfUnstructuredDataSet1Display.OpacityTransferFunction = 'PiecewiseFunction'
gradientOfUnstructuredDataSet1Display.DataAxesGrid = 'GridAxesRepresentation'
gradientOfUnstructuredDataSet1Display.SelectionCellLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.SelectionPointLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.PolarAxes = 'PolarAxesRepresentation'
gradientOfUnstructuredDataSet1Display.ScalarOpacityUnitDistance = 0.04225672412170924
gradientOfUnstructuredDataSet1Display.ScalarOpacityFunction = vorticityPWF
gradientOfUnstructuredDataSet1Display.Slice = 31

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
gradientOfUnstructuredDataSet1Display.ScaleTransferFunction.Points = [-6.663405101660004, 0.0, 0.5, 0.0, 6.663405101732195, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
gradientOfUnstructuredDataSet1Display.OpacityTransferFunction.Points = [-6.663405101660004, 0.0, 0.5, 0.0, 6.663405101732195, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
gradientOfUnstructuredDataSet1Display.DataAxesGrid.XTitleFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.YTitleFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.ZTitleFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.XLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.YLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
gradientOfUnstructuredDataSet1Display.PolarAxes.PolarAxisTitleFontFile = ''
gradientOfUnstructuredDataSet1Display.PolarAxes.PolarAxisLabelFontFile = ''
gradientOfUnstructuredDataSet1Display.PolarAxes.LastRadialAxisTextFontFile = ''
gradientOfUnstructuredDataSet1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(calculator1)
# ----------------------------------------------------------------
