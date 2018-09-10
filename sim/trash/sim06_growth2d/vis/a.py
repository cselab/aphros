# state file generated using paraview version 5.5.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1000, 500]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [1.0, 0.5, 0.125]
renderView1.StereoType = 0
renderView1.CameraPosition = [1.0, 0.5, -4.221666218300808]
renderView1.CameraFocalPoint = [1.0, 0.5, 0.125]
renderView1.CameraParallelScale = 0.5156040363942931
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

# create a new 'XDMF Reader'
vf_ = XDMFReader(FileNames=['/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_0.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_1.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_2.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_3.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_4.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_5.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_6.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_7.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_8.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_9.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_10.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_11.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_12.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_13.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_14.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_15.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_16.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_17.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_18.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_19.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_20.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_21.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_22.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_23.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_24.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_25.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_26.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_27.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_28.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_29.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_30.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_31.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_32.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_33.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_34.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_35.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_36.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_37.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_38.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_39.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_40.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_41.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_42.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_43.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_44.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_45.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_46.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_47.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_48.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_49.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_50.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_51.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_52.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_53.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_54.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_55.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_56.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_57.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_58.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_59.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_60.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_61.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_62.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_63.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_64.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_65.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_66.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_67.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_68.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_69.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_70.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_71.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_72.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_73.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_74.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_75.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_76.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_77.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_78.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_79.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_80.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_81.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_82.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_83.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_84.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_85.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_86.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_87.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_88.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_89.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_90.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_91.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_92.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_93.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_94.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_95.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_96.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_97.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_98.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_99.xmf', '/home/kpetr/s/falcon_home/mfer/ch/sim/sim06_growth2d/vf_100.xmf'])
vf_.CellArrayStatus = ['vf']
vf_.GridStatus = ['Grid_6914']

# create a new 'Outline'
outline1 = Outline(Input=vf_)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from vf_
vf_Display = Show(vf_, renderView1)

# get color transfer function/color map for 'vf'
vfLUT = GetColorTransferFunction('vf')
vfLUT.RGBPoints = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
vfLUT.ColorSpace = 'RGB'
vfLUT.NanColor = [1.0, 0.0, 0.0]
vfLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
vf_Display.Representation = 'Surface'
vf_Display.ColorArrayName = ['CELLS', 'vf']
vf_Display.LookupTable = vfLUT
vf_Display.OSPRayScaleFunction = 'PiecewiseFunction'
vf_Display.SelectOrientationVectors = 'None'
vf_Display.ScaleFactor = 0.2
vf_Display.SelectScaleArray = 'None'
vf_Display.GlyphType = 'Arrow'
vf_Display.GlyphTableIndexArray = 'None'
vf_Display.GaussianRadius = 0.01
vf_Display.SetScaleArray = [None, '']
vf_Display.ScaleTransferFunction = 'PiecewiseFunction'
vf_Display.OpacityArray = [None, '']
vf_Display.OpacityTransferFunction = 'PiecewiseFunction'
vf_Display.DataAxesGrid = 'GridAxesRepresentation'
vf_Display.SelectionCellLabelFontFile = ''
vf_Display.SelectionPointLabelFontFile = ''
vf_Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
vf_Display.DataAxesGrid.XTitleFontFile = ''
vf_Display.DataAxesGrid.YTitleFontFile = ''
vf_Display.DataAxesGrid.ZTitleFontFile = ''
vf_Display.DataAxesGrid.XLabelFontFile = ''
vf_Display.DataAxesGrid.YLabelFontFile = ''
vf_Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
vf_Display.PolarAxes.PolarAxisTitleFontFile = ''
vf_Display.PolarAxes.PolarAxisLabelFontFile = ''
vf_Display.PolarAxes.LastRadialAxisTextFontFile = ''
vf_Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from outline1
outline1Display = Show(outline1, renderView1)

# trace defaults for the display properties.
outline1Display.Representation = 'Surface'
outline1Display.ColorArrayName = [None, '']
outline1Display.LineWidth = 2.0
outline1Display.OSPRayScaleFunction = 'PiecewiseFunction'
outline1Display.SelectOrientationVectors = 'None'
outline1Display.ScaleFactor = 0.2
outline1Display.SelectScaleArray = 'None'
outline1Display.GlyphType = 'Arrow'
outline1Display.GlyphTableIndexArray = 'None'
outline1Display.GaussianRadius = 0.01
outline1Display.SetScaleArray = [None, '']
outline1Display.ScaleTransferFunction = 'PiecewiseFunction'
outline1Display.OpacityArray = [None, '']
outline1Display.OpacityTransferFunction = 'PiecewiseFunction'
outline1Display.DataAxesGrid = 'GridAxesRepresentation'
outline1Display.SelectionCellLabelFontFile = ''
outline1Display.SelectionPointLabelFontFile = ''
outline1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
outline1Display.DataAxesGrid.XTitleFontFile = ''
outline1Display.DataAxesGrid.YTitleFontFile = ''
outline1Display.DataAxesGrid.ZTitleFontFile = ''
outline1Display.DataAxesGrid.XLabelFontFile = ''
outline1Display.DataAxesGrid.YLabelFontFile = ''
outline1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
outline1Display.PolarAxes.PolarAxisTitleFontFile = ''
outline1Display.PolarAxes.PolarAxisLabelFontFile = ''
outline1Display.PolarAxes.LastRadialAxisTextFontFile = ''
outline1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'vf'
vfPWF = GetOpacityTransferFunction('vf')
vfPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(vf_)
# ----------------------------------------------------------------