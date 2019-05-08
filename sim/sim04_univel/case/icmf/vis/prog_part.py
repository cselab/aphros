import numpy as np

# execfile("/home/kpetr/s/electrochem/conf/icmf2019/slides/figloc/dim3_pltr005b000/filtc.py")

N = vtk.numpy_interface.dataset_adapter.VTKNoneArray()

A = inputs[0]

AC = A.CellData
AP = A.PointData

q = 'c'

sel = None
for D in [AC, AP]:
    if q in D.keys():
        sel = np.unique(D[q])

sel = [a for a in sel if (a % 123) % 7 == 0]

if q in AC.keys():
    cc = AC[q]
    ss = np.array([c if c in sel else -1 for c in cc])
    output.ShallowCopy(A.VTKObject)
    output.CellData.append(ss , "sc")

if q in AP.keys():
    cc = AP[q]
    ss = np.array([c if c in sel else -1 for c in cc])
    output.ShallowCopy(A.VTKObject)
    output.PointData.append(ss , "sp")
