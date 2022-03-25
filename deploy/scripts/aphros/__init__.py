from aphros.par import *
import aphros.io
from aphros.vtk import ReadVtkPoly

from aphros.io import *

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
