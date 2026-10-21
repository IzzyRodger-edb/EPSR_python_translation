"""Module data of ``module bragginten_mod`` (fortran_src/ato_manipulation_rev4.f, lines 20-29).

Allocatable arrays start as ``None``; when allocated they are 1-based arrays from
``fortran.arrays`` (index 0 unused).
"""

import sys

from . import register


def reset():
    """(Re)create the module variables of ``module bragginten_mod``."""
    global hihkl, kihkl, lihkl, ramp, iamp, intenhkl, intenqunique, qhkl, qunique
    global nhkl, ntotalqunique, nqunique, iqvalues, ibvalues
    # integer, dimension(:), allocatable :: hihkl,kihkl,lihkl
    hihkl = None
    kihkl = None
    lihkl = None
    # real, dimension(:,:), allocatable :: ramp,iamp,intenhkl,intenqunique
    # Real and imaginary parts of scattering amplitude for each atom type and Bragg reflection
    ramp = None
    iamp = None
    intenhkl = None
    intenqunique = None
    # real, dimension(:), allocatable :: qhkl,qunique      !Q values of Bragg peaks
    qhkl = None
    qunique = None
    # integer :: nhkl,ntotalqunique
    nhkl = 0
    ntotalqunique = 0
    # integer, dimension(:), allocatable :: nqunique,iqvalues
    nqunique = None
    iqvalues = None
    # integer, dimension(:,:), allocatable :: ibvalues
    ibvalues = None


register(sys.modules[__name__])
reset()
