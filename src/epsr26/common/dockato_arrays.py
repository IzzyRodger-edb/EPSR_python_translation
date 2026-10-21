"""Fortran module data of ``module dockato_arrays`` (fortran_src/dockato.f lines 1-20).

Hand-written (TRANSLATION_GUIDE 5.2).  The module has no routines; it is used by
``dockato`` and ``bond_molecule`` (programs/dockato.py) to hold the first molecule of the
docking .ato file.

Allocatable arrays start as ``None``; uninitialised module scalars are 0 (gfortran puts
module variables in zero-initialised static storage).
"""

import sys

from ..fortran.fmath import f32
from . import register


def reset():
    """(Re)create the module variables of dockato_arrays."""
    global natos, ngrots, natomax, ngrotmax, nrotmax, ndihmax, nadd
    global nbondss, ndihedralangs, iabonds, labaxis, labrotdihs
    global xadds, yadds, zadds, amassts, phixs, phiys, phizs
    global xas, yas, zas, tbonds, dihedralangs, overlapfactor
    global nrott1, nrott2, labrott1, labrott2
    global alabels, rotgrps
    natos = 0
    ngrots = 0
    natomax = 0
    ngrotmax = 0
    nrotmax = 0
    ndihmax = 0
    nadd = 0

    nbondss = None  # integer(:)
    ndihedralangs = None  # integer(:)
    iabonds = None  # integer(:,:)
    labaxis = None  # integer(:,:)
    labrotdihs = None  # integer(:,:,:)

    xadds = f32(0.0)
    yadds = f32(0.0)
    zadds = f32(0.0)
    amassts = f32(0.0)
    phixs = f32(0.0)
    phiys = f32(0.0)
    phizs = f32(0.0)
    xas = None  # real(:)
    yas = None  # real(:)
    zas = None  # real(:)
    tbonds = None  # real(:,:)
    dihedralangs = None  # real(:,:)
    overlapfactor = f32(0.0)  # Used to control the extent to which overlaps are allowed. 0.0 means all overlaps are allowed.

    # temporary rotational groups
    nrott1 = 0
    nrott2 = 0
    labrott1 = None  # integer(:)
    labrott2 = None  # integer(:)

    alabels = None  # character(len=3)(:)
    rotgrps = None  # character(len=3)(:)


register(sys.modules[__name__])
reset()
