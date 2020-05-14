from aphros.par import *
import aphros.plot
import aphros.io
import aphros.stream
from aphros.vtk import ReadVtkPoly

ReadPlain = aphros.io.ReadPlain
InitBasicFigure = aphros.plot.InitBasicFigure
SaveBasicFigure = aphros.plot.SaveBasicFigure
PlotFieldCoolwarm = aphros.plot.PlotFieldCoolwarm
stream = aphros.stream.stream

import aphros.confgen
Geometry = aphros.confgen.Geometry
Bc = aphros.confgen.Bc
Var = aphros.confgen.Var
