"""Fortran module data of ``module sdfcube_mod`` (fortran_src/sdfcube.f, lines 1-13).

Hand-written (TRANSLATION_GUIDE 5.2).  The module has no routines; the program and its
routines (``run_sdfcube``, ``initcube``, ``outcube``, ``sdfcube``) are in
``programs/sdfcube.py``.  (``common/sdfcube.py`` is the include file sdfcube.inc.)
"""

import sys

from . import register


def reset():
    """(Re)create the module variables of sdfcube_mod."""
    global natosave, lenfname_cube, xasave, yasave, zasave, atsymbolsave, rhocube
    global fname_cube, nm, phim, cthetm, chim, xcent, ycent, zcent
    global performtest, outputaveragemol
    natosave = 0
    lenfname_cube = 0
    xasave = None  # real, dimension(:), allocatable
    yasave = None
    zasave = None
    atsymbolsave = None  # character(len=3), dimension(:), allocatable
    rhocube = None  # real, dimension(:,:,:), allocatable
    fname_cube = ''  # character(len=256)
    nm = None  # integer, dimension(:), allocatable: first or second molecule
    phim = None  # real, dimension(:), allocatable
    cthetm = None
    chim = None
    xcent = None
    ycent = None
    zcent = None
    # Logical variables to perform a test calculation and to output average or single
    # molecules for plotting purposes
    performtest = False
    outputaveragemol = False


register(sys.modules[__name__])
reset()
