"""Fortran module data of ``module lngamma_mod`` (fortran_src/lngamma_mod.f90).

Hand-written (TRANSLATION_GUIDE 5.2).  The routines (``fctrln``, ``lngamma``) are in
``library/lngamma_mod.py``.
"""

import sys

from . import register


def reset():
    """(Re)create the module variables of lngamma_mod."""
    global mgamma, gamln
    mgamma = 0
    gamln = None  # real(kind=8), dimension(:), allocatable


register(sys.modules[__name__])
reset()
