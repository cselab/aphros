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
renderView1.ViewSize = [1502, 918]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.5, 0.5, 0.03125]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.533094612090724, 0.46718514573490644, 2.7659675244629605]
renderView1.CameraFocalPoint = [0.533094612090724, 0.46718514573490644, 0.03125]
renderView1.CameraParallelScale = 0.6626320255546222
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

# create a new 'XML Structured Grid Reader'
p_0 = XMLStructuredGridReader(FileName=['/home/kpetr/s/ch_partstr/univel/ch/p_0000.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0001.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0002.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0003.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0004.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0005.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0006.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0007.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0008.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0009.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0010.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0011.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0012.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0013.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0014.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0015.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0016.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0017.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0018.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0019.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0020.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0021.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0022.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0023.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0024.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0025.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0026.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0027.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0028.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0029.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0030.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0031.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0032.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0033.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0034.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0035.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0036.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0037.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0038.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0039.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0040.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0041.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0042.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0043.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0044.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0045.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0046.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0047.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0048.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0049.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0050.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0051.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0052.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0053.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0054.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0055.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0056.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0057.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0058.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0059.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0060.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0061.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0062.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0063.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0064.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0065.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0066.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0067.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0068.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0069.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0070.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0071.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0072.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0073.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0074.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0075.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0076.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0077.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0078.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0079.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0080.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0081.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0082.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0083.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0084.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0085.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0086.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0087.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0088.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0089.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0090.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0091.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0092.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0093.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0094.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0095.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0096.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0097.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0098.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0099.vts', '/home/kpetr/s/ch_partstr/univel/ch/p_0100.vts'])
p_0.CellArrayStatus = ['vx', 'vy', 'vz', 'p', 'vf', 'cl']

# create a new 'Slice'
slice1 = Slice(Input=p_0)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5, 0.5, 0.03125]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# create a new 'CSV Reader'
partit_0 = CSVReader(FileName=['/home/kpetr/s/ch_partstr/univel/ch/partit_0000.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0001.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0002.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0003.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0004.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0005.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0006.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0007.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0008.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0009.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0010.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0011.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0012.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0013.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0014.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0015.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0016.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0017.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0018.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0019.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0020.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0021.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0022.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0023.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0024.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0025.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0026.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0027.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0028.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0029.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0030.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0031.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0032.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0033.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0034.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0035.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0036.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0037.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0038.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0039.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0040.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0041.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0042.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0043.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0044.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0045.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0046.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0047.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0048.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0049.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0050.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0051.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0052.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0053.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0054.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0055.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0056.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0057.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0058.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0059.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0060.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0061.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0062.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0063.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0064.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0065.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0066.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0067.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0068.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0069.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0070.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0071.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0072.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0073.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0074.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0075.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0076.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0077.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0078.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0079.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0080.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0081.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0082.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0083.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0084.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0085.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0086.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0087.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0088.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0089.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0090.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0091.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0092.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0093.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0094.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0095.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0096.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0097.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0098.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0099.csv', '/home/kpetr/s/ch_partstr/univel/ch/partit_0100.csv'])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=partit_0)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Calculator'
calculator1 = Calculator(Input=tableToPoints1)
calculator1.ResultArrayName = 'cc'
calculator1.Function = 'sin(c*1234567)'

# create a new 'Legacy VTK Reader'
s_0 = LegacyVTKReader(FileNames=['/home/kpetr/s/ch_partstr/univel/ch/s_0000.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0001.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0002.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0003.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0004.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0005.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0006.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0007.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0008.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0009.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0010.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0011.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0012.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0013.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0014.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0015.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0016.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0017.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0018.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0019.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0020.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0021.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0022.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0023.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0024.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0025.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0026.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0027.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0028.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0029.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0030.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0031.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0032.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0033.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0034.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0035.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0036.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0037.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0038.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0039.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0040.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0041.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0042.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0043.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0044.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0045.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0046.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0047.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0048.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0049.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0050.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0051.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0052.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0053.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0054.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0055.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0056.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0057.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0058.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0059.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0060.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0061.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0062.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0063.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0064.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0065.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0066.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0067.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0068.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0069.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0070.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0071.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0072.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0073.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0074.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0075.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0076.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0077.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0078.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0079.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0080.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0081.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0082.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0083.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0084.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0085.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0086.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0087.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0088.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0089.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0090.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0091.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0092.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0093.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0094.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0095.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0096.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0097.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0098.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0099.vtk', '/home/kpetr/s/ch_partstr/univel/ch/s_0100.vtk'])

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from p_0
p_0Display = Show(p_0, renderView1)

# trace defaults for the display properties.
p_0Display.Representation = 'Wireframe'
p_0Display.ColorArrayName = ['POINTS', '']
p_0Display.Opacity = 0.53
p_0Display.OSPRayScaleFunction = 'PiecewiseFunction'
p_0Display.SelectOrientationVectors = 'None'
p_0Display.ScaleFactor = 0.1
p_0Display.SelectScaleArray = 'None'
p_0Display.GlyphType = 'Arrow'
p_0Display.GlyphTableIndexArray = 'None'
p_0Display.GaussianRadius = 0.005
p_0Display.SetScaleArray = [None, '']
p_0Display.ScaleTransferFunction = 'PiecewiseFunction'
p_0Display.OpacityArray = [None, '']
p_0Display.OpacityTransferFunction = 'PiecewiseFunction'
p_0Display.DataAxesGrid = 'GridAxesRepresentation'
p_0Display.SelectionCellLabelFontFile = ''
p_0Display.SelectionPointLabelFontFile = ''
p_0Display.PolarAxes = 'PolarAxesRepresentation'
p_0Display.ScalarOpacityUnitDistance = 0.2229420780051279

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
p_0Display.DataAxesGrid.XTitleFontFile = ''
p_0Display.DataAxesGrid.YTitleFontFile = ''
p_0Display.DataAxesGrid.ZTitleFontFile = ''
p_0Display.DataAxesGrid.XLabelFontFile = ''
p_0Display.DataAxesGrid.YLabelFontFile = ''
p_0Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
p_0Display.PolarAxes.PolarAxisTitleFontFile = ''
p_0Display.PolarAxes.PolarAxisLabelFontFile = ''
p_0Display.PolarAxes.LastRadialAxisTextFontFile = ''
p_0Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from slice1
slice1Display = Show(slice1, renderView1)

# get color transfer function/color map for 'vf'
vfLUT = GetColorTransferFunction('vf')
vfLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.5000000000000001, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]
vfLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'vf']
slice1Display.LookupTable = vfLUT
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.005
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.SelectionCellLabelFontFile = ''
slice1Display.SelectionPointLabelFontFile = ''
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
slice1Display.DataAxesGrid.XTitleFontFile = ''
slice1Display.DataAxesGrid.YTitleFontFile = ''
slice1Display.DataAxesGrid.ZTitleFontFile = ''
slice1Display.DataAxesGrid.XLabelFontFile = ''
slice1Display.DataAxesGrid.YLabelFontFile = ''
slice1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
slice1Display.PolarAxes.PolarAxisTitleFontFile = ''
slice1Display.PolarAxes.PolarAxisLabelFontFile = ''
slice1Display.PolarAxes.LastRadialAxisTextFontFile = ''
slice1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from s_0
s_0Display = Show(s_0, renderView1)

# trace defaults for the display properties.
s_0Display.Representation = 'Wireframe'
s_0Display.ColorArrayName = [None, '']
s_0Display.LineWidth = 3.0
s_0Display.OSPRayScaleFunction = 'PiecewiseFunction'
s_0Display.SelectOrientationVectors = 'None'
s_0Display.ScaleFactor = 0.041070280224084856
s_0Display.SelectScaleArray = 'None'
s_0Display.GlyphType = 'Arrow'
s_0Display.GlyphTableIndexArray = 'None'
s_0Display.GaussianRadius = 0.0020535140112042426
s_0Display.SetScaleArray = [None, '']
s_0Display.ScaleTransferFunction = 'PiecewiseFunction'
s_0Display.OpacityArray = [None, '']
s_0Display.OpacityTransferFunction = 'PiecewiseFunction'
s_0Display.DataAxesGrid = 'GridAxesRepresentation'
s_0Display.SelectionCellLabelFontFile = ''
s_0Display.SelectionPointLabelFontFile = ''
s_0Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
s_0Display.DataAxesGrid.XTitleFontFile = ''
s_0Display.DataAxesGrid.YTitleFontFile = ''
s_0Display.DataAxesGrid.ZTitleFontFile = ''
s_0Display.DataAxesGrid.XLabelFontFile = ''
s_0Display.DataAxesGrid.YLabelFontFile = ''
s_0Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
s_0Display.PolarAxes.PolarAxisTitleFontFile = ''
s_0Display.PolarAxes.PolarAxisLabelFontFile = ''
s_0Display.PolarAxes.LastRadialAxisTextFontFile = ''
s_0Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from calculator1
calculator1Display = Show(calculator1, renderView1)

# get color transfer function/color map for 'cc'
ccLUT = GetColorTransferFunction('cc')
ccLUT.RGBPoints = [-0.9951680157308399, 0.0, 0.0, 1.0, -0.6655392952313994, 0.0, 0.0, 1.0, -0.663553580047668, 1.0, 0.0, 1.0, -0.33591057473195907, 1.0, 0.0, 1.0, -0.3339248595482276, 0.0, 1.0, 1.0, -0.002310423865055644, 0.0, 1.0, 1.0, -0.0003247086813240596, 0.0, 1.0, 0.0, 0.3273182966343847, 0.0, 1.0, 0.0, 0.3293040118181164, 1.0, 1.0, 0.0, 0.6569470171338251, 1.0, 1.0, 0.0, 0.6589327323175566, 1.0, 0.0, 0.0, 0.9905471680007286, 1.0, 0.0, 0.0]
ccLUT.ColorSpace = 'HSV'
ccLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['POINTS', 'cc']
calculator1Display.LookupTable = ccLUT
calculator1Display.PointSize = 5.0
calculator1Display.LineWidth = 3.0
calculator1Display.RenderPointsAsSpheres = 1
calculator1Display.OSPRayScaleArray = 'c'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 0.04113186073913205
calculator1Display.SelectScaleArray = 'None'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'None'
calculator1Display.GaussianRadius = 0.0020565930369566025
calculator1Display.SetScaleArray = ['POINTS', 'c']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'c']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.SelectionCellLabelFontFile = ''
calculator1Display.SelectionPointLabelFontFile = ''
calculator1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [3004.0, 0.0, 0.5, 0.0, 5005.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [3004.0, 0.0, 0.5, 0.0, 5005.0, 1.0, 0.5, 0.0]

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

# setup the color legend parameters for each legend in this view

# get color legend/bar for vfLUT in view renderView1
vfLUTColorBar = GetScalarBar(vfLUT, renderView1)
vfLUTColorBar.Title = 'vf'
vfLUTColorBar.ComponentTitle = ''
vfLUTColorBar.TitleFontFile = ''
vfLUTColorBar.LabelFontFile = ''

# set color bar visibility
vfLUTColorBar.Visibility = 1

# get color legend/bar for ccLUT in view renderView1
ccLUTColorBar = GetScalarBar(ccLUT, renderView1)
ccLUTColorBar.WindowLocation = 'UpperRightCorner'
ccLUTColorBar.Title = 'cc'
ccLUTColorBar.ComponentTitle = ''
ccLUTColorBar.TitleFontFile = ''
ccLUTColorBar.LabelFontFile = ''

# set color bar visibility
ccLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'vf'
vfPWF = GetOpacityTransferFunction('vf')
vfPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'cc'
ccPWF = GetOpacityTransferFunction('cc')
ccPWF.Points = [-0.9951680157308399, 0.0, 0.5, 0.0, 0.9905471680007286, 1.0, 0.5, 0.0]
ccPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(calculator1)
# ----------------------------------------------------------------