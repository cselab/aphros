from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

renderView1 = GetRenderView()
renderView1.InteractionMode = '2D'

if 'base' not in globals():
    base = 'sides4'

nx032 = LegacyVTKReader(
    registrationName='nx032',
    FileNames=[
        base + '_nx032/s_0000.vtk', base + '_nx032/s_0010.vtk',
        base + '_nx032/s_0020.vtk', base + '_nx032/s_0030.vtk',
        base + '_nx032/s_0040.vtk', base + '_nx032/s_0050.vtk',
        base + '_nx032/s_0060.vtk', base + '_nx032/s_0070.vtk',
        base + '_nx032/s_0080.vtk'
    ])

nx064 = LegacyVTKReader(
    registrationName='nx064',
    FileNames=[
        base + '_nx064/s_0000.vtk', base + '_nx064/s_0010.vtk',
        base + '_nx064/s_0020.vtk', base + '_nx064/s_0030.vtk',
        base + '_nx064/s_0040.vtk', base + '_nx064/s_0050.vtk',
        base + '_nx064/s_0060.vtk', base + '_nx064/s_0070.vtk',
        base + '_nx064/s_0080.vtk'
    ])

nx128 = LegacyVTKReader(
    registrationName='nx128',
    FileNames=[
        base + '_nx128/s_0000.vtk', base + '_nx128/s_0010.vtk',
        base + '_nx128/s_0020.vtk', base + '_nx128/s_0030.vtk',
        base + '_nx128/s_0040.vtk', base + '_nx128/s_0050.vtk',
        base + '_nx128/s_0060.vtk', base + '_nx128/s_0070.vtk',
        base + '_nx128/s_0080.vtk'
    ])

nx256 = LegacyVTKReader(
    registrationName='nx256',
    FileNames=[
        base + '_nx256/s_0000.vtk', base + '_nx256/s_0010.vtk',
        base + '_nx256/s_0020.vtk', base + '_nx256/s_0030.vtk',
        base + '_nx256/s_0040.vtk', base + '_nx256/s_0050.vtk',
        base + '_nx256/s_0060.vtk', base + '_nx256/s_0070.vtk',
        base + '_nx256/s_0080.vtk'
    ])

nx512 = LegacyVTKReader(
    registrationName='nx512',
    FileNames=[
        base + '_nx512/s_0000.vtk', base + '_nx512/s_0010.vtk',
        base + '_nx512/s_0020.vtk', base + '_nx512/s_0030.vtk',
        base + '_nx512/s_0040.vtk', base + '_nx512/s_0050.vtk',
        base + '_nx512/s_0060.vtk', base + '_nx512/s_0070.vtk',
        base + '_nx512/s_0080.vtk'
    ])

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

nx032Display = Show(nx032, renderView1)
nx032Display.Representation = 'Wireframe'
nx032Display.ColorArrayName = ['POINTS', '']

nx064Display = Show(nx064, renderView1)
nx064Display.Representation = 'Wireframe'
nx064Display.ColorArrayName = ['POINTS', '']

nx128Display = Show(nx128, renderView1)
nx128Display.Representation = 'Wireframe'
nx128Display.ColorArrayName = ['POINTS', '']

nx256Display = Show(nx256, renderView1)
nx256Display.Representation = 'Wireframe'
nx256Display.ColorArrayName = ['POINTS', '']

nx512Display = Show(nx512, renderView1)
nx512Display.Representation = 'Wireframe'
nx512Display.ColorArrayName = ['POINTS', '']

