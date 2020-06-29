import sys
import re
import math
import os
import numpy as np
import paraview.simple as para

def Log(s):
    sys.stderr.write(str(s) + '\n')

def GetStep(path):
    return re.findall('[^_]*_([0-9]*)\.*.', os.path.basename(path))[0]

def GetSteps(paths):
    return list(map(GetStep, paths))

def ReplaceFilename(paths, pattern, keep_dir=True):
    """
    Replaces filename by pattern with step index substitution.
    paths: `list(str)`
        Paths.
    pattern: `str`
        Pattern containing a single `{}` to be replaced by step index.

    Example:
    >>> ReplaceFilename(["dir/vf_0001.xmf"], "sm_{}.vtk")
    >>> ["dir/sm_0001.vtk"]
    """
    r = []
    for f in paths:
        dirname = os.path.dirname(f)
        basename = os.path.basename(f)
        step = GetStep(f)
        if keep_dir:
            r.append(os.path.join(dirname, pattern.format(step)))
        else:
            r.append(pattern.format(step))
    return r

def SubstituteStep(*args, **kwargs):
    Log("Warning: SubstituteStep() is deprecated. Renamed to ReplaceFilename()")
    return ReplaceFilename(*args, **kwargs)


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

def SaveAnimation(steps, renderView, sources_ft, timearrays, pattern="a_{:}.png", force=False):
    for index,step in enumerate(steps):
        outfile = pattern.format(step)
        if os.path.isfile(outfile) and not force:
            Log("skip existing {:}".format(outfile))
            continue
        SetTimeStep(index, sources_ft, timearrays)
        Log("{:}/{:}: {:}".format(index + 1, len(steps), outfile))
        para.SaveScreenshot(outfile, renderView)
