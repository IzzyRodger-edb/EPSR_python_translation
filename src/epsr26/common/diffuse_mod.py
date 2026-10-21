"""Module data of ``module diffuse_mod`` (fortran_src/ato_manipulation_rev4.f, lines 38-111).

The module routine ``initialise_diffuse_values`` is in ``library/ato_manipulation_rev4.py``.
Allocatable arrays start as ``None``.
"""

import sys

from ..fortran.fmath import f32
from . import register


def reset():
    """(Re)create the module variables of ``module diffuse_mod``."""
    global multdiffa, multdiffb, multdiffc, ndiffuse, nqdiff, ldiffmax, ldiffstep, mdiffstep
    global l1ex, n1ex, ndiffcoeff, diffcoeffc, diffcoeffs
    global diffcoefftc, diffcoeffts, diffcoeffoldc, diffcoeffolds, qdiff, qdiffstep
    global diffusedefined, diffusename
    # integer :: multdiffa,multdiffb,multdiffc,ndiffuse,nqdiff,ldiffmax,ldiffstep,mdiffstep
    # !No. of unit cells along each of h,k,l
    multdiffa = 0
    multdiffb = 0
    multdiffc = 0
    ndiffuse = 0
    nqdiff = 0
    ldiffmax = 0
    ldiffstep = 0
    mdiffstep = 0
    # integer, dimension(:), allocatable :: l1ex,n1ex,ndiffcoeff ! Sets the pattern along each axis
    l1ex = None
    n1ex = None
    ndiffcoeff = None
    # real, dimension(:,:,:), allocatable :: diffcoeffc,diffcoeffs
    diffcoeffc = None
    diffcoeffs = None
    # real, dimension(:,:), allocatable :: diffcoefftc,diffcoeffts,diffcoeffoldc,diffcoeffolds
    diffcoefftc = None
    diffcoeffts = None
    diffcoeffoldc = None
    diffcoeffolds = None
    # real, dimension(:), allocatable :: qdiff
    qdiff = None
    # real :: qdiffstep
    qdiffstep = f32(0.0)
    # logical :: diffusedefined
    diffusedefined = False
    # character(len=3) :: diffusename
    diffusename = ''


register(sys.modules[__name__])
reset()
