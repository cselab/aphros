# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

traj = GetActiveSource()

tableToPoints1 = TableToPoints(Input=traj)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

view = FindView("RenderView1")

tableToPoints1Display = Show(tableToPoints1, view)

glyph1 = Glyph(Input=tableToPoints1, GlyphType='Sphere')
glyph1.GlyphType = 'Sphere'
glyph1.ScaleArray = ['POINTS', 'r']
glyph1.ScaleFactor = 1.0
glyph1.GlyphMode = 'All Points'

glyph1Display = Show(glyph1, view)

glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
glyph1Display.ColorArrayName = [None, '']
