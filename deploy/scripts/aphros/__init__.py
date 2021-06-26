from aphros.par import *
import aphros.plot
import aphros.io
import aphros.stream
from aphros.vtk import ReadVtkPoly

from aphros.io import *
InitBasicFigure = aphros.plot.InitBasicFigure
SaveBasicFigure = aphros.plot.SaveBasicFigure
PlotFieldCoolwarm = aphros.plot.PlotFieldCoolwarm

InitFigure = aphros.plot.InitFigure
SaveFigure = aphros.plot.SaveFigure
ApplyParams = aphros.plot.ApplyParams
PlotField = aphros.plot.PlotField
PlotFieldText = aphros.plot.PlotFieldText

ReplaceFilename = aphros.plot.ReplaceFilename
GetStep = aphros.plot.GetStep
GetSteps = aphros.plot.GetSteps
stream = aphros.stream.stream

import aphros.confgen
Geometry = aphros.confgen.Geometry
BoundaryConditions = aphros.confgen.BoundaryConditions
Bc = aphros.confgen.BoundaryConditions
Config = aphros.confgen.Config
ReadConfig = aphros.confgen.ReadConfig
Var = aphros.confgen.Config
Domain = aphros.confgen.Domain
AdjustedDomain = aphros.confgen.AdjustedDomain
PartitionDomain = aphros.confgen.PartitionDomain
GenerateJobConfig = aphros.confgen.GenerateJobConfig
AdjustDomainToProcessors = aphros.confgen.AdjustDomainToProcessors
CheckDomain = aphros.confgen.CheckDomain
from aphros.confgen import Parameters

from aphros.test import TestBase
