#!/usr/bin/env python3

import argparse
import math
import os
import sys
import copy
import subprocess
try:
    import numpy as np
except ImportError:
    pass

def NormalizeType(v):
    if isinstance(v, (float, np.float32, np.float64)):
        v = float(v)
    elif isinstance(v, (int, np.int32, np.int64)):
        v = int(v)
    elif isinstance(v, str):
        v = str(v)
    elif isinstance(v, list):
        v = list(map(float, v))
    else:
        assert False, "unknown type of value '{:}'".format(str(v))
    return v


class Config:
    def __init__(self, var=dict()):
        for k, v in var.items():
            self[k] = v

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def Generate(self):
        t = []
        for k, v in vars(self).items():
            v = NormalizeType(v)
            if isinstance(v, float):
                t.append("set double {:} {:.20g}".format(k, v))
            elif isinstance(v, int):
                t.append("set int {:} {:}".format(k, v))
            elif isinstance(v, list):
                v = ' '.join(["{:.20g}".format(a) for a in v])
                t.append("set vect {:} {:}".format(k, v))
            elif isinstance(v, str):
                if '\n' in v:
                    v = '"' + v + '"'
                t.append("set string {:} {:}".format(k, v))
            else:
                assert False, "unknown type of value '{:}'".format(str(v))
        return '\n'.join(t)

    def GenerateFile(self, path):
        with open(path, 'w') as f:
            f.write(self.Generate())


def VectToStr(v):
    return ' '.join(map(str, v))


class Domain:
    """
    Description of the domain size, the number of blocks and processors.
    All arguments are optional but must provide enough
    information to define th domain size and mesh size.

    lx,ly,lz: `float`
        Domain length.
    nx,ny,nz: `int`
        Number of cells in each direction.
    h: `float`
        Cell length.
    extent: `float`
        Maximum domain length over all directions.
    bsx,bsy,bsz: `int`
        Block size in each direction.
    bx,by,bz: `int`
        Number of blocks in each direction.
    px,py,pz: `int`
        Number of processors in each direction.
    nproc: `int`
        Number of processors.
    """
    def __init__(self,
                 lx=None,
                 ly=None,
                 lz=None,
                 nx=None,
                 ny=None,
                 nz=None,
                 h=None,
                 extent=None,
                 bsx=None,
                 bsy=None,
                 bsz=None,
                 bs=32,
                 bx=None,
                 by=None,
                 bz=None,
                 nproc=None):
        conf = locals()

        def tryline(line):
            nonlocal conf
            try:
                exec(line, None, conf)
            except:
                pass

        def trylines(lines):
            for line in lines:
                tryline(line)

        # TODO: detect conflicting parameters

        trylines([
            "extent = max(lx, ly, lz)",
            "nmax = max(nx, ny, nz)",
            "bsx = bs if bsx is None else bsx",
            "bsy = bs if bsy is None else bsy",
            "bsz = bs if bsz is None else bsz",
            "bs = bsx",
            "h = lz / nz",
            "h = lx / nx",
            "h = ly / ny",
            "h = extent / nmax",
            "lx = float(nx * h)",
            "ly = float(ny * h)",
            "lz = float(nz * h)",
            "nx = int(lx / h + 0.5)",
            "ny = int(ly / h + 0.5)",
            "nz = int(lz / h + 0.5)",
            "bx = (nx + bsx - 1) // bsx",
            "by = (ny + bsy - 1) // bsy",
            "bz = (nz + bsz - 1) // bsz",
            "nx = bx * bsx",
            "ny = by * bsy",
            "nz = bz * bsz",
            "lx = float(nx * h)",
            "ly = float(ny * h)",
            "lz = float(nz * h)",
            "extent = max(lx, ly, lz)",
            "nmax = max(nx, ny, nz)",
            "px = 1",
            "py = 1",
            "pz = 1",
        ])

        for k, v in conf.items():
            if k in ['self']:
                continue
            setattr(self, k, v)

    def __str__(self):
        return "Domain: {}".format(str(vars(self)))

    def compare(self, other):
        def isclose(a, b):
            return np.all(np.isclose(a, b))

        diff = dict()
        v = vars(self)
        vo = vars(other)
        for k in v:
            if not isclose(v[k], vo[k]):
                diff[k] = (v[k], vo[k])
        return diff

    def printdiff(diff):
        for k, v in diff.items():
            print("{:>8}  {:<7} {:<7}".format(k, *v))

    def GetMeshConfig(self):
        return """\
set int px {px}
set int py {py}
set int pz {pz}

set int bx {bx}
set int by {by}
set int bz {bz}

set int bsx {bsx}
set int bsy {bsy}
set int bsz {bsz}
""".format(**vars(self))

    def GenerateMeshConfig(self, path="mesh.conf"):
        with open(path, 'w') as f:
            f.write(self.GetMeshConfig())


def PartitionDomain(domain):
    """
    Partitions domain to subdomains for given number of processors
    minimizing the communication area.
    Returns new domain with updated px,py,pz and divided bx,by,bz.
    """
    CheckDomain(domain)
    dim = 3
    b_all = [domain.bx, domain.by, domain.bz]
    bs = [domain.bsx, domain.bsy, domain.bsz]

    def prod(v):
        res = 1
        for a in v:
            res *= a
        return res

    def divisible(p, d):
        return all(p[i] % d[i] == 0 for i in range(dim))

    def div(p, d):
        return [p[i] // d[i] for i in range(dim)]

    def quality(p):
        nonlocal bs, b_all
        b = div(b_all, p)
        return -sum([prod(b) // b[i] * prod(bs) // bs[i] for i in range(dim)])

    nproc = domain.nproc
    divisors = [i for i in range(1, nproc + 1) if nproc % i == 0]
    best_p = None
    for px in divisors:
        for py in divisors:
            pz = nproc // (px * py)
            p = [px, py, pz]
            if prod(p) != nproc:
                continue
            if not divisible(b_all, p):
                continue
            if best_p is None or quality(p) > quality(best_p):
                best_p = p
    assert best_p, "No partition found for {:}".format(domain)
    newdomain = copy.deepcopy(domain)
    p = best_p
    newdomain.px = p[0]
    newdomain.py = p[1]
    newdomain.pz = p[2]
    b = div(b_all, p)
    newdomain.bx = b[0]
    newdomain.by = b[1]
    newdomain.bz = b[2]
    return newdomain


def AdjustedDomain(verbose=True,
                   decrease_nproc=True,
                   increase_size=False,
                   nproc_factor=0.9,
                   size_factor=1.1,
                   **kwargs):
    domain = Domain(**kwargs)
    newdomain = AdjustDomainToProcessors(
        domain,
        decrease_nproc=decrease_nproc,
        increase_size=increase_size,
        nproc_factor=nproc_factor,
        size_factor=size_factor,
    )
    assert newdomain, "Compatible domain not found for \n{:}".format(domain)
    diff = domain.compare(newdomain)
    if diff and verbose:
        print("Domain is adjusted. Changes:")
        Domain.printdiff(diff)
    CheckDomain(newdomain)
    return newdomain


def CheckDomain(domain, fatal=True):
    """
    Reports statistics and performs quality checks of the domain:
        - number of blocks is divisible by `nproc`

    Returns list of failed checks or aborts if fatal=True.
    """
    checks = [
        "(domain.bx * domain.by * domain.bz) % domain.nproc == 0",
    ]
    failed = []
    for check in checks:
        if not eval(check):
            failed.append(check)
    assert not failed, \
            "While checking domain\n{:}\nthe following checks failed\n{:}".format(
        domain, failed)
    return failed


def AdjustDomainToProcessors(domain,
                             decrease_nproc=True,
                             increase_size=False,
                             nproc_factor=0.9,
                             size_factor=1.1,
                             verbose=True):
    """
    Adjusts the domain size and the number of processors
    to make the number of blocks divisible by `domain.nproc`.

    domain: `Domain`
        Domain to start from.
    decrease_nproc: `bool`
        Allow less processors, down to `domain.nproc * nproc_factor`.
    increase_size: `bool`
        Allow larger domain, multiply the number
        of blocks up to `size_factor`.
    nproc_factor: `float`
        Factor for minimal allowed number of processors.

    Returns:
    newdomain: `Domain`
        Adjusted domain or None if not found.
    """
    b_range = lambda b: range(b, int(b * size_factor + 0.5) + 1)
    nproc_range = lambda n: range(n, int(n * nproc_factor + 0.5) - 1, -1)
    for bz in b_range(domain.bz):
        for by in b_range(domain.by):
            for bx in b_range(domain.bx):
                for nproc in nproc_range(domain.nproc):
                    if (bx * by * bz) % nproc == 0:
                        d = Domain(bx=bx,
                                   by=by,
                                   bz=bz,
                                   h=domain.h,
                                   bsx=domain.bsx,
                                   bsy=domain.bsy,
                                   bsz=domain.bsz,
                                   bs=domain.bs,
                                   nproc=nproc)
                        return d
    return None


class Geometry:
    def __init__(self):
        self.lines = []

    def __Append(self, line):
        self.lines.append(line)

    def __Prefix(self, intersect=False, invert=False):
        s = ''
        if intersect:
            s += '&'
        if invert:
            s += '-'
        return s

    def Box(self, center, halfsize, rotation_z=0, **kwargs):
        s = self.__Prefix(**kwargs)
        s += "box {:}   {:}   {:}".format(VectToStr(center),
                                          VectToStr(halfsize), rotation_z)
        self.__Append(s)
        return self

    def RoundBox(self, center, halfsize, roundradius=0, **kwargs):
        s = self.__Prefix(**kwargs)
        s += "roundbox {:}   {:}   {:}".format(VectToStr(center),
                                               VectToStr(halfsize),
                                               roundradius)
        self.__Append(s)
        return self

    def Sphere(self, center, radii, **kwargs):
        s = self.__Prefix(**kwargs)
        s += "sphere {:}   {:}".format(VectToStr(center), VectToStr(radii))
        self.__Append(s)
        return self

    def Cylinder(self, center, normal, radius, normalrange, **kwargs):
        s = self.__Prefix(**kwargs)
        s += "cylinder {:}   {:}   {:}   {:}".format(VectToStr(center),
                                                     VectToStr(normal), radius,
                                                     VectToStr(normalrange))
        self.__Append(s)
        return self

    def Polygon(self, origin, normal, right, normalrange, scale, polygon,
                **kwargs):
        '''
        polygon: `list(list(float))`, polygon as list of 2D vertices
            If first and last vertices do not coincide,
            appended by the first vertex to make a loop
        '''
        assert all(len(p) == 2 for p in polygon), \
                "expected 2D points, got '{:}'".format(polygon)
        assert len(polygon) >= 3, \
                "expected polygon of at least 3 vertices, got '{:}'".format(polygon)
        if polygon[0] != polygon[-1]:
            polygon.append(polygon[0])
        coords = [x for p in polygon for x in p]
        s = self.__Prefix(**kwargs)
        s += "polygon {:}   {:}   {:}   {:}   {:}   {:}".format(
            VectToStr(origin), VectToStr(normal), VectToStr(right),
            VectToStr(normalrange), scale, VectToStr(coords))
        self.__Append(s)
        return self

    def Polygon2(self, origin, right, scale, polygon, **kwargs):
        '''
        polygon: `list(list(float))`, polygon as list of 2D vertices
            If first and last vertices do not coincide,
            appended by the first vertex to make a loop
        '''
        assert all(len(p) == 2 for p in polygon), \
                "expected 2D points, got '{:}'".format(polygon)
        assert len(polygon) >= 3, \
                "expected polygon of at least 3 vertices, got '{:}'".format(polygon)
        if polygon[0] != polygon[-1]:
            polygon.append(polygon[0])
        coords = [x for p in polygon for x in p]
        s = self.__Prefix(**kwargs)
        s += "polygon2 {:}   {:}   {:}   {:}".format(VectToStr(origin),
                                                     VectToStr(right), scale,
                                                     VectToStr(coords))
        self.__Append(s)
        return self

    def Ruled(self, origin, normal, right, normalrange, scale0, scale1,
              polygon0, polygon1, **kwargs):
        '''
        polygon0,polygon1: `list(list(float))`
            Polygon as lists of 2D vertices used on the opposite
            sides of the ruled surface scaled by `scale0` and `scale1`
            respectively.
        '''
        for polygon in [polygon0, polygon1]:
            assert all(len(p) == 2 for p in polygon), \
                    "expected 2D points, got '{:}'".format(polygon)
            assert len(polygon) >= 3, \
                    "expected polygons of at least 3 vertices, got '{:}'".format(polygon)
            if polygon[0] != polygon[-1]:
                polygon.append(polygon[0])
        coords = [x for p in polygon0 + polygon1 for x in p]
        s = self.__Prefix(**kwargs)
        s += "ruled {:}   {:}   {:}   {:}   {:} {:}   {:}".format(
            VectToStr(origin), VectToStr(normal), VectToStr(right),
            VectToStr(normalrange), scale0, scale1, VectToStr(coords))
        self.__Append(s)
        return self

    def Generate(self):
        return '\n'.join(self.lines)

    def GenerateFile(self, path):
        with open(path, 'w') as f:
            f.write(self.Generate())


class BoundaryConditions:
    def __init__(self):
        self.lines = []

    def __Indent(self, text):
        lines = text.split('\n')
        return '\n'.join(["  " + line for line in lines])

    def __Append(self, line, geom):
        line = "{:} {{\n{:}\n}}".format(line, self.__Indent(geom.Generate()))
        self.lines.append(line)

    def Wall(self, geom, velocity=[0, 0, 0], extra=""):
        s = "wall {:}".format(VectToStr(velocity)) + extra
        self.__Append(s, geom)

    def WallRotation(self, geom, center, omega):
        s = "wall_rotation {:} {:}".format(VectToStr(center), VectToStr(omega))
        self.__Append(s, geom)

    def WallRotationMagn(self, geom, center, omega):
        s = "wall_rotation_magn {:} {:}".format(VectToStr(center),
                                                VectToStr(omega))
        self.__Append(s, geom)

    def SlipWall(self, geom, extra=""):
        s = "slipwall" + extra
        self.__Append(s, geom)

    def Inlet(self, geom, velocity, extra=""):
        s = "inlet {:}{}".format(VectToStr(velocity), extra)
        self.__Append(s, geom)

    def InletFlux(self, geom, velocity, index):
        s = "inletflux {:} {:}".format(VectToStr(velocity), index)
        self.__Append(s, geom)

    def InletPressure(self, geom, pressure, extra=""):
        s = "inletpressure {:}{}".format(pressure, extra)
        self.__Append(s, geom)

    def InletRotation(self, geom, veln, center, omega, extra=""):
        s = "inlet_rotation {:}   {:}   {:}{}".format(veln, VectToStr(center),
                                                      VectToStr(omega), extra)
        self.__Append(s, geom)

    def Outlet(self, geom, extra=""):
        s = "outlet" + extra
        self.__Append(s, geom)

    def OutletPressure(self, geom, pressure, extra=""):
        s = "outletpressure {:}{}".format(pressure, extra)
        self.__Append(s, geom)

    def Symm(self, geom):
        s = "symm"
        self.__Append(s, geom)

    def Custom(self, geom, desc):
        self.__Append(desc, geom)

    def Generate(self):
        return '\n'.join(self.lines)

    def GenerateFile(self, path):
        with open(path, 'w') as f:
            f.write(self.Generate())


class Parameters:
    """
    Base class for parameters of a config generator.
    Values of parameters are accessible as attributes:
    >>> par = Parameters()
    >>> par.extent = 1
    """
    def __init__(self, path=None):
        for k in dir(self):
            if not k.startswith('_') and k not in ['exec']:
                setattr(self, k, getattr(self, k))
        if path:
            self.execfile(path)

    def execfile(self, path):
        """
        Executes a Python script and saves the local variables as attributes.
        path: `str`
            Path to script.
        """
        with open(path) as f:
            d = vars(self)
            exec(f.read(), None, d)
            for k, v in d.items():
                setattr(self, k, v)
        return self


def ReadConfig(fpath):
    """
    Returns Config from configuration file (commands "set ...").
    fpath: path to configuration file
    """
    code = subprocess.check_output(['ap.conf2py', fpath])
    d = dict()
    exec(code, None, d)
    return Config(d)

def GenerateJobConfig(nproc, time_limit_minutes, basedir="."):
    with open(os.path.join(basedir, "np"), 'w') as f:
        f.write(str(nproc))
    with open(os.path.join(basedir, "tl"), 'w') as f:
        f.write(str(time_limit_minutes))
