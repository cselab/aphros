import numpy
import math
import sys

av = sys.argv
# length of polygon in normal direction
LN = float(av[0])
# length of polygon in tangential direction
LT = float(av[1])
# displacement in normal direction
DN = float(av[2])
# offset in particle index (+np for next string)
OI = int(av[3])
# offset in particle index (+np for next string)
TYPE = av[4] # 'poly' or 'line'

if TYPE == "poly":
   poly = True
elif TYPE == "line":
   poly = False
else:
   assert False

def Cross(a, b):
   return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
        ]

def P(x, y, z):
   return numpy.array((x, y, z))

def Norm(p):
   return sum(p * p) ** 0.5

# o: origin (central particle)
# pa,pb: points (adjacent particle)
def GetPlane(o, pa, pb):
   global Norm
   global LT, LN, DN
   en = (pa + pb - 2 * o)
   en /= Norm(en)
   et = pb - pa
   et /= Norm(et)
   rr = [[-1,-1], [1,-1], [1,1], [-1,1]]
   pp = [o + et * t * LT * 0.5 + en * (n * LN * 0.5 + DN) for t,n in rr]
   return pp

# o: origin (central particle)
# pa,pb: points (adjacent particle)
def GetLine(o, pa, pb):
   global Norm
   global LT, LN, DN
   en = (pa + pb - 2 * o)
   en /= Norm(en)
   et = pb - pa
   et /= Norm(et)
   pp = [o + en * (LN * 0.5 + DN), o + en * (-LN * 0.5 + DN)]
   return pp


a = inputs[0]
ap = a.Points
ad = a.PointData
ii, = numpy.where(ad['sp'] >= 0)

ns = 2  # number of strings per cell
np = len(ii) // ns

assert len(ii) % ns == 0

nph = np // 2
nph += OI

pp = vtk.vtkPoints()

so = ap[ii[nph]]
sop = ap[ii[nph + 1]]
som = ap[ii[nph - 1]]

so = P(*so)
som = P(*som)
sop = P(*sop)

pl = GetPlane(so, som, sop) if poly else GetLine(so, som, sop)
n = len(pl)

for i in range(n):
   pp.InsertPoint(i, *pl[i])

o = self.GetPolyDataOutput()
o.SetPoints(pp)

ll = vtk.vtkPolygon() if poly else vtk.vtkPolyLine()

ll.GetPointIds().SetNumberOfIds(n)
for i in range(n):
   ll.GetPointIds().SetId(i, i)

o.Allocate(1, 1)

o.InsertNextCell(ll.GetCellType(), ll.GetPointIds())
