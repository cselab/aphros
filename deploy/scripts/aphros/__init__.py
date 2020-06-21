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
BoundaryConditions = aphros.confgen.BoundaryConditions
Bc = aphros.confgen.BoundaryConditions
Config = aphros.confgen.Config
Var = aphros.confgen.Config
Domain = aphros.confgen.Domain
AdjustedDomain = aphros.confgen.AdjustedDomain
AdjustDomainToProcessors = aphros.confgen.AdjustDomainToProcessors
CheckDomain = aphros.confgen.CheckDomain

from aphros.test import TestBase
