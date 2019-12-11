import sys
import os

interactive = "__file__" not in globals()
if interactive:
    assert "vtkpython" in sys.executable

filename = ""
if interactive:
    import PyQt5.QtGui
    import PyQt5.QtWidgets
    print(dir(PyQt5.QtGui))
    filename = PyQt5.QtWidgets.QFileDialog.getOpenFileNames(None, "Open input files")

me = os.path.dirname(sys.argv[0])
d = os.path.dirname(os.path.realpath(me))
d = "/tmp"
with open(os.path.join(d, 'o'), 'w') as f:
    def W(s):
        f.write(str(s) + '\n')
    W(["interactive", interactive])
    W(me)
    W(os.path.realpath(me))
    W(sys.executable)
    W(sys.argv)
    W(os.getcwd())
    W(["filename", filename])
