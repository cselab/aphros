from glob import glob
import sys
import re
import math
import os
import numpy as np
import paraview.simple as para

def Log(s):
    sys.stderr.write(str(s) + '\n')

av = sys.argv

if len(av) < 2:
    sys.stderr.write('''usage: {:} [datafiles]'''.format(
            os.path.basename(av[0])))
    exit(1)

files = av[1:]
basenames = list(map(os.path.basename, files))
dirnames = list(map(os.path.dirname, files))
steps = [re.findall("_([0-9]*)", f)[0] for f in basenames]

def GetArgFiles():
    '''
    Returns list of files passed as arguments.
    '''
    return files

def GetArgSteps():
    '''
    Returns step indices from filenames in arguments.
    '''
    return steps

def GetArgDir():
    '''
    Returns directory containing the first file in arguments.
    '''
    return dirnames[0]

def FindFiles(pattern=None, basedir=GetArgDir(), check=True):
    '''
    Returns list of files with names given by pattern.
    '''
    ff = GetArgFiles() if pattern is None else \
            [os.path.join(basedir, pattern.format(s)) for s in GetArgSteps()]

    if check:
        for f in ff:
            assert os.path.isfile(f), "file not found: '{:}'".format(f)
    return ff

def ApplyForceTime(sources):
    '''
    Applies ForceTime filter to sources.
    Returns new sources and arrays with original time values.
    '''
    timearrays = [np.array(s.TimestepValues) for s in sources]
    sources_ft = [para.ForceTime(s) for s in sources]
    return sources_ft, timearrays

def SetTimeStep(index, sources_ft, timearrays):
    '''
    Sets sources to time step given by index.
    '''
    for s,tt in zip(sources_ft, timearrays):
        s.ForcedTime = tt[index]
        s.UpdatePipeline()

def GetBox(o):
    '''
    Returns bounding box of object o.
    '''
    o.UpdatePipeline()
    di = o.GetDataInformation()
    lim = di.DataInformation.GetBounds()
    lim0 = np.array(lim[::2])
    lim1 = np.array(lim[1::2])
    return lim0, lim1

def SaveAnimation(renderView, sources_ft, timearrays, pattern="a_{:}.png"):
    steps = GetArgSteps()
    for index,step in enumerate(steps):
        fout = pattern.format(step)
        if os.path.isfile(fout):
            Log("skip existing {:}".format(fout))
            continue
        SetTimeStep(index, sources_ft, timearrays)
        Log("{:}/{:}: {:}".format(index + 1, len(steps), fout))
        para.SaveScreenshot(fout, renderView)
