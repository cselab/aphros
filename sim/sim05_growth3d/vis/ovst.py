# state file generated using paraview version 5.5.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Light'
light1 = CreateLight()
light1.Intensity = 0.8
light1.Position = [4.0, 0.5, -8.0]
light1.FocalPoint = [4.0, 0.5, -6.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# create light
# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1536, 1186]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [4.0, 0.498749986290932, 0.524999976158142]
renderView1.StereoType = 0
renderView1.CameraPosition = [7.05083065037415, 1.99328509007997, -15.0767923414364]
renderView1.CameraFocalPoint = [4.0, 0.498749986290932, 0.524999976158142]
renderView1.CameraViewUp = [-0.0010369827893544996, 0.9954619246947966, 0.09515503743694166]
renderView1.CameraParallelScale = 3.58547823281106
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.196078431372549, 0.196078431372549, 0.196078431372549]
renderView1.EnableOSPRay = 1
renderView1.Shadows = 1
renderView1.AmbientSamples = 10
renderView1.SamplesPerPixel = 10
renderView1.AdditionalLights = light1
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
platey0 = Box()
platey0.XLength = 8.0
platey0.YLength = 0.05
platey0.Center = [4.0, -0.025, 0.5]

# create a new 'XDMF Reader'
vf_ = XDMFReader(FileNames=['/home/kpetr/s/euler_sim06/a1tdblal/vf_0.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_1.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_2.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_3.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_4.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_5.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_6.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_7.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_8.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_9.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_10.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_11.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_12.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_13.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_14.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_15.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_16.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_17.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_18.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_19.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_20.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_21.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_22.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_23.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_24.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_25.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_26.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_27.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_28.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_29.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_30.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_31.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_32.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_33.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_34.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_35.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_36.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_37.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_38.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_39.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_40.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_41.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_42.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_43.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_44.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_45.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_46.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_47.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_48.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_49.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_50.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_51.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_52.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_53.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_54.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_55.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_56.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_57.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_58.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_59.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_60.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_61.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_62.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_63.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_64.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_65.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_66.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_67.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_68.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_69.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_70.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_71.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_72.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_73.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_74.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_75.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_76.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_77.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_78.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_79.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_80.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_81.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_82.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_83.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_84.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_85.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_86.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_87.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_88.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_89.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_90.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_91.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_92.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_93.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_94.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_95.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_96.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_97.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_98.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_99.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_100.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_101.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_102.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_103.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_104.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_105.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_106.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_107.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_108.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_109.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_110.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_111.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_112.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_113.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_114.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_115.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_116.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_117.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_118.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_119.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_120.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_121.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_122.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_123.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_124.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_125.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_126.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_127.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_128.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_129.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_130.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_131.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_132.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_133.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_134.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_135.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_136.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_137.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_138.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_139.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_140.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_141.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_142.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_143.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_144.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_145.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_146.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_147.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_148.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_149.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_150.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_151.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_152.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_153.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_154.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_155.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_156.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_157.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_158.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_159.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_160.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_161.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_162.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_163.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_164.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_165.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_166.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_167.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_168.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_169.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_170.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_171.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_172.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_173.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_174.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_175.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_176.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_177.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_178.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_179.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_180.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_181.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_182.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_183.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_184.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_185.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_186.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_187.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_188.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_189.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_190.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_191.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_192.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_193.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_194.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_195.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_196.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_197.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_198.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_199.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_200.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_201.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_202.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_203.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_204.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_205.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_206.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_207.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_208.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_209.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_210.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_211.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_212.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_213.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_214.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_215.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_216.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_217.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_218.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_219.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_220.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_221.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_222.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_223.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_224.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_225.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_226.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_227.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_228.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_229.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_230.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_231.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_232.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_233.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_234.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_235.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_236.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_237.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_238.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_239.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_240.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_241.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_242.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_243.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_244.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_245.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_246.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_247.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_248.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_249.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_250.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_251.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_252.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_253.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_254.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_255.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_256.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_257.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_258.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_259.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_260.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_261.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_262.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_263.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_264.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_265.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_266.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_267.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_268.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_269.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_270.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_271.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_272.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_273.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_274.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_275.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_276.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_277.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_278.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_279.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_280.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_281.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_282.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_283.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_284.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_285.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_286.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_287.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_288.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_289.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_290.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_291.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_292.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_293.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_294.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_295.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_296.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_297.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_298.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_299.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_300.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_301.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_302.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_303.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_304.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_305.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_306.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_307.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_308.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_309.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_310.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_311.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_312.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_313.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_314.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_315.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_316.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_317.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_318.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_319.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_320.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_321.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_322.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_323.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_324.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_325.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_326.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_327.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_328.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_329.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_330.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_331.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_332.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_333.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_334.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_335.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_336.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_337.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_338.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_339.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_340.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_341.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_342.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_343.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_344.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_345.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_346.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_347.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_348.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_349.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_350.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_351.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_352.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_353.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_354.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_355.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_356.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_357.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_358.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_359.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_360.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_361.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_362.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_363.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_364.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_365.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_366.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_367.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_368.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_369.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_370.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_371.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_372.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_373.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_374.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_375.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_376.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_377.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_378.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_379.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_380.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_381.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_382.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_383.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_384.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_385.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_386.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_387.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_388.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_389.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_390.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_391.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_392.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_393.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_394.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_395.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_396.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_397.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_398.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_399.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_400.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_401.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_402.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_403.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_404.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_405.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_406.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_407.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_408.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_409.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_410.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_411.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_412.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_413.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_414.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_415.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_416.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_417.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_418.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_419.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_420.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_421.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_422.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_423.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_424.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_425.xmf', '/home/kpetr/s/euler_sim06/a1tdblal/vf_426.xmf'])
vf_.CellArrayStatus = ['vf']
vf_.GridStatus = ['Grid_11552']

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=vf_)

# create a new 'Box'
electrode0 = Box()
electrode0.XLength = 0.05
electrode0.YLength = 0.4
electrode0.ZLength = 0.05
electrode0.Center = [0.525, -0.1975, 0.075]

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=vf_)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingDimensions = [512, 64, 67]
resampleToImage1.SamplingBounds = [0.0, 8.0, 0.0, 1.0, -0.05, 1.0]

# create a new 'Box'
electrode1 = Box()
electrode1.XLength = 0.05
electrode1.YLength = 0.4
electrode1.ZLength = 0.05
electrode1.Center = [0.525, 1.195, 0.075]

# create a new 'Append Geometry'
electrodes = AppendGeometry(Input=[electrode0, electrode1])

# create a new 'Contour'
contour1 = Contour(Input=resampleToImage1)
contour1.ContourBy = ['POINTS', 'vf']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Clip'
oxygen = Clip(Input=contour1)
oxygen.ClipType = 'Plane'
oxygen.Scalars = ['POINTS', 'vtkValidPointMask']

# init the 'Plane' selected for 'ClipType'
oxygen.ClipType.Origin = [3.28865742683411, 0.5, 0.134546503424644]
oxygen.ClipType.Normal = [0.0, -1.0, 0.0]

# create a new 'Legacy VTK Reader'
s = LegacyVTKReader(FileNames=['/home/kpetr/s/euler_sim06/a1tdblal/s.0.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.1.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.2.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.3.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.4.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.5.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.6.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.7.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.8.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.9.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.10.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.11.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.12.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.13.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.14.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.15.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.16.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.17.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.18.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.19.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.20.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.21.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.22.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.23.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.24.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.25.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.26.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.27.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.28.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.29.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.30.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.31.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.32.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.33.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.34.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.35.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.36.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.37.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.38.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.39.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.40.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.41.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.42.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.43.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.44.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.45.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.46.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.47.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.48.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.49.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.50.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.51.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.52.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.53.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.54.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.55.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.56.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.57.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.58.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.59.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.60.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.61.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.62.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.63.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.64.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.65.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.66.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.67.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.68.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.69.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.70.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.71.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.72.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.73.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.74.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.75.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.76.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.77.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.78.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.79.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.80.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.81.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.82.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.83.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.84.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.85.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.86.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.87.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.88.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.89.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.90.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.91.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.92.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.93.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.94.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.95.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.96.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.97.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.98.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.99.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.100.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.101.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.102.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.103.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.104.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.105.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.106.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.107.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.108.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.109.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.110.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.111.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.112.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.113.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.114.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.115.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.116.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.117.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.118.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.119.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.120.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.121.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.122.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.123.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.124.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.125.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.126.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.127.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.128.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.129.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.130.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.131.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.132.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.133.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.134.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.135.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.136.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.137.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.138.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.139.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.140.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.141.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.142.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.143.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.144.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.145.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.146.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.147.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.148.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.149.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.150.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.151.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.152.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.153.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.154.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.155.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.156.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.157.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.158.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.159.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.160.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.161.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.162.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.163.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.164.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.165.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.166.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.167.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.168.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.169.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.170.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.171.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.172.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.173.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.174.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.175.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.176.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.177.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.178.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.179.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.180.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.181.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.182.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.183.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.184.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.185.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.186.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.187.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.188.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.189.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.190.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.191.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.192.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.193.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.194.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.195.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.196.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.197.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.198.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.199.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.200.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.201.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.202.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.203.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.204.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.205.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.206.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.207.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.208.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.209.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.210.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.211.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.212.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.213.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.214.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.215.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.216.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.217.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.218.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.219.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.220.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.221.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.222.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.223.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.224.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.225.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.226.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.227.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.228.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.229.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.230.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.231.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.232.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.233.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.234.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.235.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.236.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.237.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.238.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.239.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.240.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.241.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.242.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.243.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.244.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.245.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.246.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.247.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.248.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.249.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.250.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.251.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.252.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.253.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.254.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.255.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.256.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.257.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.258.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.259.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.260.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.261.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.262.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.263.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.264.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.265.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.266.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.267.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.268.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.269.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.270.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.271.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.272.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.273.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.274.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.275.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.276.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.277.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.278.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.279.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.280.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.281.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.282.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.283.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.284.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.285.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.286.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.287.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.288.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.289.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.290.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.291.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.292.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.293.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.294.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.295.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.296.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.297.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.298.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.299.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.300.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.301.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.302.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.303.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.304.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.305.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.306.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.307.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.308.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.309.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.310.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.311.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.312.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.313.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.314.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.315.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.316.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.317.vtk', '/home/kpetr/s/euler_sim06/a1tdblal/s.318.vtk'])

# create a new 'Box'
platez1 = Box()
platez1.XLength = 8.0
platez1.YLength = 1.1
platez1.ZLength = 0.05
platez1.Center = [4.0, 0.5, 1.025]

# create a new 'Clip'
hydrogen = Clip(Input=contour1)
hydrogen.ClipType = 'Plane'
hydrogen.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
hydrogen.ClipType.Origin = [3.28865742683411, 0.5, 0.134546503424644]
hydrogen.ClipType.Normal = [0.0, 1.0, 0.0]

# create a new 'Box'
platez0 = Box()
platez0.XLength = 8.0
platez0.YLength = 1.1
platez0.ZLength = 0.05
platez0.Center = [4.0, 0.5, -0.025]

# create a new 'Box'
platey1 = Box()
platey1.XLength = 8.0
platey1.YLength = 0.05
platey1.Center = [4.0, 1.025, 0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from platey0
platey0Display = Show(platey0, renderView1)

# trace defaults for the display properties.
platey0Display.Representation = 'Surface'
platey0Display.ColorArrayName = [None, '']
platey0Display.DiffuseColor = [0.454901960784314, 0.454901960784314, 0.454901960784314]
platey0Display.OSPRayScaleArray = 'Normals'
platey0Display.OSPRayScaleFunction = 'PiecewiseFunction'
platey0Display.SelectOrientationVectors = 'None'
platey0Display.ScaleFactor = 0.8
platey0Display.SelectScaleArray = 'None'
platey0Display.GlyphType = 'Arrow'
platey0Display.GlyphTableIndexArray = 'None'
platey0Display.GaussianRadius = 0.00899999991059303
platey0Display.SetScaleArray = ['POINTS', 'Normals']
platey0Display.ScaleTransferFunction = 'PiecewiseFunction'
platey0Display.OpacityArray = ['POINTS', 'Normals']
platey0Display.OpacityTransferFunction = 'PiecewiseFunction'
platey0Display.DataAxesGrid = 'GridAxesRepresentation'
platey0Display.SelectionCellLabelFontFile = ''
platey0Display.SelectionPointLabelFontFile = ''
platey0Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
platey0Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
platey0Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
platey0Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.XTitleFontFile = ''
platey0Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.YTitleFontFile = ''
platey0Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.ZTitleFontFile = ''
platey0Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.XLabelFontFile = ''
platey0Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.YLabelFontFile = ''
platey0Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
platey0Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
platey0Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
platey0Display.PolarAxes.PolarAxisTitleFontFile = ''
platey0Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
platey0Display.PolarAxes.PolarAxisLabelFontFile = ''
platey0Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
platey0Display.PolarAxes.LastRadialAxisTextFontFile = ''
platey0Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
platey0Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from electrodes
electrodesDisplay = Show(electrodes, renderView1)

# trace defaults for the display properties.
electrodesDisplay.Representation = 'Surface'
electrodesDisplay.ColorArrayName = [None, '']
electrodesDisplay.DiffuseColor = [0.835294117647059, 0.556862745098039, 0.0]
electrodesDisplay.Diffuse = 0.63
electrodesDisplay.OSPRayScaleArray = 'Normals'
electrodesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
electrodesDisplay.SelectOrientationVectors = 'None'
electrodesDisplay.ScaleFactor = 0.8
electrodesDisplay.SelectScaleArray = 'None'
electrodesDisplay.GlyphType = 'Arrow'
electrodesDisplay.GlyphTableIndexArray = 'None'
electrodesDisplay.GaussianRadius = 0.00899999991059303
electrodesDisplay.SetScaleArray = ['POINTS', 'Normals']
electrodesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
electrodesDisplay.OpacityArray = ['POINTS', 'Normals']
electrodesDisplay.OpacityTransferFunction = 'PiecewiseFunction'
electrodesDisplay.DataAxesGrid = 'GridAxesRepresentation'
electrodesDisplay.SelectionCellLabelFontFile = ''
electrodesDisplay.SelectionPointLabelFontFile = ''
electrodesDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
electrodesDisplay.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
electrodesDisplay.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
electrodesDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.XTitleFontFile = ''
electrodesDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.YTitleFontFile = ''
electrodesDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.ZTitleFontFile = ''
electrodesDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.XLabelFontFile = ''
electrodesDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.YLabelFontFile = ''
electrodesDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
electrodesDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
electrodesDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
electrodesDisplay.PolarAxes.PolarAxisTitleFontFile = ''
electrodesDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
electrodesDisplay.PolarAxes.PolarAxisLabelFontFile = ''
electrodesDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
electrodesDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
electrodesDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
electrodesDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from platey1
platey1Display = Show(platey1, renderView1)

# trace defaults for the display properties.
platey1Display.Representation = 'Surface'
platey1Display.AmbientColor = [0.0, 0.0, 0.0]
platey1Display.ColorArrayName = [None, '']
platey1Display.DiffuseColor = [0.454901960784314, 0.454901960784314, 0.454901960784314]
platey1Display.OSPRayScaleArray = 'Normals'
platey1Display.OSPRayScaleFunction = 'PiecewiseFunction'
platey1Display.SelectOrientationVectors = 'None'
platey1Display.ScaleFactor = 0.8
platey1Display.SelectScaleArray = 'None'
platey1Display.GlyphType = 'Arrow'
platey1Display.GlyphTableIndexArray = 'None'
platey1Display.GaussianRadius = 0.00899999991059303
platey1Display.SetScaleArray = ['POINTS', 'Normals']
platey1Display.ScaleTransferFunction = 'PiecewiseFunction'
platey1Display.OpacityArray = ['POINTS', 'Normals']
platey1Display.OpacityTransferFunction = 'PiecewiseFunction'
platey1Display.DataAxesGrid = 'GridAxesRepresentation'
platey1Display.SelectionCellLabelFontFile = ''
platey1Display.SelectionPointLabelFontFile = ''
platey1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
platey1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
platey1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
platey1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.XTitleFontFile = ''
platey1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.YTitleFontFile = ''
platey1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.ZTitleFontFile = ''
platey1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.XLabelFontFile = ''
platey1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.YLabelFontFile = ''
platey1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
platey1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
platey1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
platey1Display.PolarAxes.PolarAxisTitleFontFile = ''
platey1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
platey1Display.PolarAxes.PolarAxisLabelFontFile = ''
platey1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
platey1Display.PolarAxes.LastRadialAxisTextFontFile = ''
platey1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
platey1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from platez1
platez1Display = Show(platez1, renderView1)

# trace defaults for the display properties.
platez1Display.Representation = 'Surface'
platez1Display.AmbientColor = [0.0, 0.0, 0.0]
platez1Display.ColorArrayName = [None, '']
platez1Display.DiffuseColor = [0.454901960784314, 0.454901960784314, 0.454901960784314]
platez1Display.OSPRayScaleArray = 'Normals'
platez1Display.OSPRayScaleFunction = 'PiecewiseFunction'
platez1Display.SelectOrientationVectors = 'None'
platez1Display.ScaleFactor = 0.8
platez1Display.SelectScaleArray = 'None'
platez1Display.GlyphType = 'Arrow'
platez1Display.GlyphTableIndexArray = 'None'
platez1Display.GaussianRadius = 0.00899999991059303
platez1Display.SetScaleArray = ['POINTS', 'Normals']
platez1Display.ScaleTransferFunction = 'PiecewiseFunction'
platez1Display.OpacityArray = ['POINTS', 'Normals']
platez1Display.OpacityTransferFunction = 'PiecewiseFunction'
platez1Display.DataAxesGrid = 'GridAxesRepresentation'
platez1Display.SelectionCellLabelFontFile = ''
platez1Display.SelectionPointLabelFontFile = ''
platez1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
platez1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
platez1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
platez1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.XTitleFontFile = ''
platez1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.YTitleFontFile = ''
platez1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.ZTitleFontFile = ''
platez1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.XLabelFontFile = ''
platez1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.YLabelFontFile = ''
platez1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
platez1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
platez1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
platez1Display.PolarAxes.PolarAxisTitleFontFile = ''
platez1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
platez1Display.PolarAxes.PolarAxisLabelFontFile = ''
platez1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
platez1Display.PolarAxes.LastRadialAxisTextFontFile = ''
platez1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
platez1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from hydrogen
hydrogenDisplay = Show(hydrogen, renderView1)

# trace defaults for the display properties.
hydrogenDisplay.Representation = 'Surface'
hydrogenDisplay.AmbientColor = [0.0, 0.0, 0.0]
hydrogenDisplay.ColorArrayName = [None, '']
hydrogenDisplay.DiffuseColor = [0.541176470588235, 0.541176470588235, 0.541176470588235]
hydrogenDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
hydrogenDisplay.SelectOrientationVectors = 'None'
hydrogenDisplay.ScaleFactor = 0.549072283506393
hydrogenDisplay.SelectScaleArray = 'None'
hydrogenDisplay.GlyphType = 'Arrow'
hydrogenDisplay.GlyphTableIndexArray = 'None'
hydrogenDisplay.GaussianRadius = 0.0274536141753197
hydrogenDisplay.SetScaleArray = [None, '']
hydrogenDisplay.ScaleTransferFunction = 'PiecewiseFunction'
hydrogenDisplay.OpacityArray = [None, '']
hydrogenDisplay.OpacityTransferFunction = 'PiecewiseFunction'
hydrogenDisplay.DataAxesGrid = 'GridAxesRepresentation'
hydrogenDisplay.SelectionCellLabelFontFile = ''
hydrogenDisplay.SelectionPointLabelFontFile = ''
hydrogenDisplay.PolarAxes = 'PolarAxesRepresentation'
hydrogenDisplay.ScalarOpacityUnitDistance = 0.312750485185418

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
hydrogenDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.XTitleFontFile = ''
hydrogenDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.YTitleFontFile = ''
hydrogenDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.ZTitleFontFile = ''
hydrogenDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.XLabelFontFile = ''
hydrogenDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.YLabelFontFile = ''
hydrogenDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
hydrogenDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
hydrogenDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
hydrogenDisplay.PolarAxes.PolarAxisTitleFontFile = ''
hydrogenDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
hydrogenDisplay.PolarAxes.PolarAxisLabelFontFile = ''
hydrogenDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
hydrogenDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
hydrogenDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
hydrogenDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from oxygen
oxygenDisplay = Show(oxygen, renderView1)

# trace defaults for the display properties.
oxygenDisplay.Representation = 'Surface'
oxygenDisplay.AmbientColor = [0.0, 0.0, 0.0]
oxygenDisplay.ColorArrayName = [None, '']
oxygenDisplay.DiffuseColor = [0.541176470588235, 0.388235294117647, 0.392156862745098]
oxygenDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
oxygenDisplay.SelectOrientationVectors = 'None'
oxygenDisplay.ScaleFactor = 0.549072283506393
oxygenDisplay.SelectScaleArray = 'None'
oxygenDisplay.GlyphType = 'Arrow'
oxygenDisplay.GlyphTableIndexArray = 'None'
oxygenDisplay.GaussianRadius = 0.0274536141753197
oxygenDisplay.SetScaleArray = [None, '']
oxygenDisplay.ScaleTransferFunction = 'PiecewiseFunction'
oxygenDisplay.OpacityArray = [None, '']
oxygenDisplay.OpacityTransferFunction = 'PiecewiseFunction'
oxygenDisplay.DataAxesGrid = 'GridAxesRepresentation'
oxygenDisplay.SelectionCellLabelFontFile = ''
oxygenDisplay.SelectionPointLabelFontFile = ''
oxygenDisplay.PolarAxes = 'PolarAxesRepresentation'
oxygenDisplay.ScalarOpacityUnitDistance = 0.312750485185418

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
oxygenDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.XTitleFontFile = ''
oxygenDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.YTitleFontFile = ''
oxygenDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.ZTitleFontFile = ''
oxygenDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.XLabelFontFile = ''
oxygenDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.YLabelFontFile = ''
oxygenDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
oxygenDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
oxygenDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
oxygenDisplay.PolarAxes.PolarAxisTitleFontFile = ''
oxygenDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
oxygenDisplay.PolarAxes.PolarAxisLabelFontFile = ''
oxygenDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
oxygenDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
oxygenDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
oxygenDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(resampleToImage1)
# ----------------------------------------------------------------

